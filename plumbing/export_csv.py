import pandas as pd
from loguru import logger
logger.add(f"{__file__}.log", level="DEBUG")

from escalate_plumber import json_load, WF3Data, Material, EscalateCategories, Reaction, SamplerConvexHull

PlumbingData = json_load("PlumbingData.json.gz")
ReactionsValid = PlumbingData["ReactionsValid"]
GroupDict = PlumbingData["GroupDict"]
FeatureDict = PlumbingData["FeatureDict"]
WF3Entries = PlumbingData["WF3Entries"]
MaterialInventory = PlumbingData["MaterialInventory"]
ConvexHull = PlumbingData["ConvexHull"]

_useless_features = ["_feat__index", ]


def wf3data_to_record(wf3d: WF3Data, featurize=True):
    record = dict()
    for k, v in wf3d.as_dict().items():
        if k.startswith("@"):
            continue
        if isinstance(v, dict):
            cat = k
            for i_inchikey, inchikey in enumerate(sorted(v.keys())):
                molarity = v[inchikey]
                mat = Material.select_from(inchikey, MaterialInventory)
                chem_index = "{}___{}".format(cat, i_inchikey)
                record["{}___inchikey".format(chem_index)] = inchikey
                record["{}___inchi".format(chem_index)] = mat.inchi
                record["{}___chemname".format(chem_index)] = mat.name
                record["{}___molarity".format(chem_index)] = molarity
                record["{}___molarity_max".format(chem_index)] = mat.pure_molarity
                assert molarity <= mat.pure_molarity
                if featurize:
                    for featname, featval in FeatureDict[inchikey].items():
                        if featname not in _useless_features:
                            record["{}___{}".format(chem_index, featname)] = featval
        else:
            record[k] = v
    return record


if __name__ == '__main__':
    group_key = "expver-3.0_1%1%1%1"
    group_identifiers = GroupDict[group_key]

    export_hull = True

    for featurize in (False, True):
        logger.warning(f">>> Export with features: {featurize}")
        records = [wf3data_to_record(WF3Entries[i], featurize=featurize) for i in group_identifiers]

        if export_hull:
            ConvexHull.df.drop_duplicates().to_csv("csv/convexhull.csv", index=False)
            possible_inchikeys = set()
            for r in records:
                for k, v in r.items():
                    if k.endswith("___inchikey"):
                        possible_inchikeys.add(v)
            possible_inchikeys = sorted(possible_inchikeys)
            export_hull = False

        df = pd.DataFrame.from_records(records)

        # an ugly padding scheme for `1%1%1%1` reactions
        zerofill_fields = []
        valuefill_fields = []
        valuefill_table = dict()
        for col in df.columns:
            if df[col].isna().any():
                if col.endswith("___molarity"):
                    zerofill_fields.append(col)
                else:
                    valuefill_fields.append(col)
                    assert any(col.startswith(cat) for cat in EscalateCategories)
                    notna_values = sorted(set([v for v in df[col] if not pd.isna(v)]))
                    assert len(notna_values) == 1, f"cannot pad na for this column: {col}," \
                                                   f" requires one unique value, but we have {notna_values}"
                    valuefill_table[col] = notna_values[0]
        for col in zerofill_fields:
            df.loc[:, col] = df.loc[:, col].fillna(value=0.)
            logger.info(f"fill na in {col} with zeros")
        for col in valuefill_fields:
            df[col] = df.loc[:, col].fillna(value=valuefill_table[col])
            logger.info(f"fill na in {col} with {valuefill_table[col]}")
        if featurize:
            df = df.select_dtypes(["number"])
            df.to_csv("csv/{}-features.csv".format(group_key), index=False)
        else:
            df.to_csv("csv/{}.csv".format(group_key), index=False)