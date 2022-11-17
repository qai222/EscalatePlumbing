Tasks
---

## [20221117_export_reactions](20221117_export_reactions):
##### description
```
export all reactions that
a. done at HC
b. has one type of A (amine) cation
```
##### results
- [reactions__one_organic.csv](./20221117_export_reactions/reactions__one_organic.csv)
  - this contains all reactions labelled with `HC` and include only one type of `organic`
  - the fingerprint of a reaction comes from the number of unique organic/inorganic/solvent/acid chemicals
  - molarities are calculated based on the alpha vial volume

## [20221117_export_reagents](20221117_export_reagents)
##### description
```
export previously prepared reagents in a human-readable format
```
##### results
- [reaction_distinguishable](20221117_export_reagents/reaction_distinguishable)
  - a list of `csv` files, each contains the reagents of a particular chemical system (an unordered combination of chemicals)
  - rows are distinguishable by the combination of their `reaction_id` and `reagent_index`
    - for example, pure water dispensed to vial `X1` is different from the same water dispensed to vial `X2`, 
    even they come from the same bottle
- [reaction_indistinguishable](20221117_export_reagents/reaction_indistinguishable)
  - same as above except this time reagents are indistinguishable if their chemical compositions are identical
