"""
2022-11-17
Alex:
export previously prepared reagents in a human-readable format
"""

import pandas as pd

from escalate_plumber import json_load, Reaction, sort_and_group

PlumbingData = json_load("../../plumbing/PlumbingData.json.gz")
ReactionsValid = PlumbingData["ReactionsValid"]
ReactionsValid: dict[str, Reaction]

if __name__ == '__main__':
    reagent_records = []
    for rid, r in ReactionsValid.items():
        reagent_records += r.export_reagent_records(include_reaction_info=True)

    reagent_records_readable = []
    for k, g in sort_and_group(reagent_records, lambda x: x["chemsys"]).items():
        df_records = pd.DataFrame.from_records(g)
        df_records.to_csv(f"reaction_distinguishable/REAGENT@@{k}.csv", index=False)

        # do not distinguish reagent used in different reactions
        df_records_short = df_records[
            [c for c in df_records.columns if not c.startswith("reaction_") and not c.startswith("reagent_")]]
        df_records_short = df_records_short.loc[df_records_short.round(6).drop_duplicates().index]
        df_records_short.to_csv(f"reaction_indistinguishable/REAGENT@@{k}.csv", index=False)
