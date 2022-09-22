Escalate plumbing
---

A set of lightweight dataclasses for reactions defined in `ESCALATE`, specifically for
`worflow==3.x` reactions (Vapor diffusion crystallization).

### Requirement
`pip install numpy monty loguru`

### Vapor diffusion crystallization

- For each reaction, two vials are used:
    - alpha vial: a mixture of `ammonium halide`, `inorganic halide`, `solvent` and `acid`.
    - beta vial: the **antisolvent**, usually `dichloromethane`.
- Two vials are placed in a chamber at elevated temperature to allow vapor diffusion from beta vial to alpha vial. 
- The outcome is the `crystal score` from inspecting the alpha vial at the end of this reaction
  - `crystal score==1` clear solution without any solid;
  - `crystal score==2` fine powder;
  - `cyrstal score==3` small crystallites (average crystal dimension <0.1 mm)
  - `cyrstal score==4` large crystallites (average crystal dimension>0.1 mm)

For more information, see the [WF3 Protocol](https://docs.google.com/document/d/1_8FQqtNb_axzeLTTk2FwlB-SZZc4QFV3aMCrYNNcieo/). 

### Datasets

#### 1. Raw dataset 
Parameters of `918` reactions are listed in the [raw dataset](plumbing/csv/expver-3.0_1%251%251%251.csv).
For each reaction in this dataset, only one type of `ammonium halide`/`inorganic halide`/`solvent`/`acid` is used.
This csv file contains the following columns:
- `identifier`: the identifier for a reaction
- `fingerprint`: number of unique `ammonium halide`/`inorganic halide`/`solvent`/`acid`  used in this reaction
- `outcome`: the `crystal score` described above
- `organic___0___inchikey`: the `inchikey` of `ammonium halide`
- `organic___0___inchi`:    the `inchi string` of `ammonium halide`
- `organic___0___chemname`: the name of `ammonium halide`
- `organic___0___molarity`: the molarity of `ammonium halide` in alpha vial
- `inorganic___0___inchikey`: the `inchikey` of `inorganic halide`
- `inorganic___0___inchi`:    the `inchi string` of `inorganic halide`
- `inorganic___0___chemname`: the name of `inorganic halide`
- `inorganic___0___molarity`: the molarity of `inorganic halide` in alpha vial
- `solvent___0___inchikey`:   the `inchikey` of `solvent`
- `solvent___0___inchi`:      the `inchi string` of `solvent`
- `solvent___0___chemname`:   the name of `solvent`
- `solvent___0___molarity`:   the molarity of `solvent` in the alpha vial
- `acid___0___inchikey`: the `inchikey` of `acid`
- `acid___0___inchi`:    the `inchi string` of `acid`
- `acid___0___chemname`: the name of `acid`
- `acid___0___molarity`: the molarity of `acid` in the alpha vial
- `alpha_vial_volume`: the total volume (in liter) of alpha vial right at the beginning
- `beta_vial_volume`: the total volume (in liter) of beta vial right at the beginning
- `reaction_time`: reaction time
- `reaction_temperature`: reaction temperature
- `antisolvent_identity`: the `inchikey` of `antisolvent` (the only chemical in beta vial)

#### 2. Featurized dataset 
A more ML-friendly `csv` is also available [link](plumbing/csv/expver-3.0_1%251%251%251-features.csv).
This csv contains all numeric columns of the [raw dataset](plumbing/csv/expver-3.0_1%251%251%251.csv), as well as 
a list of expert-selected molecular descriptors.
All molecular descriptors are calculated using `cxcalc` from [chemaxon](https://docs.chemaxon.com/display/docs/cxcalc-calculator-functions.md).
For more information, see [this document](https://ndownloader.figstatic.com/files/35712904).
