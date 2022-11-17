"""
2022-11-17
Alex:
export all reactions that
a. done at HC
b. has one type of A (amine) cation
"""
from collections import Counter

import pandas as pd

from escalate_plumber import json_load, Reaction

PlumbingData = json_load("../../plumbing/PlumbingData.json.gz")
ReactionsValid = PlumbingData["ReactionsValid"]
WF3Entries = PlumbingData["WF3Entries"]
MaterialInventory = PlumbingData["MaterialInventory"]
ReactionsValid: dict[str, Reaction]

HC_tag = "HC"

tags = []
for r in ReactionsValid:
    tags += sorted([t for t in r.split("_") if "H" in t])
assert len(set(tags)) == 1 and tags[0] == HC_tag


def filter_reactions(tag=HC_tag, cat="organic", ncat=1) -> list[Reaction]:
    reactions = {k: v for k, v in ReactionsValid.items() if tag in k.split("_")}
    reactions: dict[str, Reaction]

    filtered = []
    for r in reactions.values():
        cats = [mat.category for mat in r.inchikey_to_material.values()]
        if Counter(cats)[cat] == ncat:
            filtered.append(r)
    return filtered


if __name__ == '__main__':
    filtered_reactions = filter_reactions(HC_tag, "organic", 1)
    records = []
    for r in filtered_reactions:
        wf3data = WF3Entries[r.identifier]
        records.append(
            wf3data.as_olympus_record(mat_inventory=MaterialInventory, feat_dict=None)
        )
    df = pd.DataFrame.from_records(records)
    human_readable_columns = [
        c for c in df.columns if "inchi" not in c and
                                 not c.endswith("molarity_max") and
                                 not c.startswith("reaction_") and
                                 "antisolvent" not in c and
                                 "volume" not in c
    ]
    df[human_readable_columns].to_csv("reactions__one_organic.csv", index=False)
