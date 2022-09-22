from __future__ import annotations

import re

import pandas as pd
from loguru import logger
from pandas._typing import FilePath

from escalate_plumber.schema import Reaction, Reagent, ReagentMaterial, Material

_precision = 5
_eps = 1 ** -_precision
_target_column = "_out_crystalscore"
_experiment_version_column = "_raw_expver"
_max_num_reagents = 10
_max_num_rm = 10

COL_TARGET = "target"
COL_NAME = "name"
COL_USELESS = "useless"
COL_FEATURE = "feature"
COL_CALCULATE = "calculate"
COL_REAGENT = "reagent"
COL_REACTION = "reaction"
COL_RAW_MMOL_INCHI_KEY = "raw_mmol_inchi_key"
COL_RAW_MOLARITY_INCHI_KEY = "raw_molarity_inchi_key"
COL_RAW_INCHIKEY_CATEGORY = "raw_inchikey_type"
COL_RAW_OTHER = "raw_other"


class ReactionParserError(Exception):
    pass


def get_subrecords(record: dict, keys: list[str]):
    reagent_record = dict()
    reaction_record = dict()
    feature_record = dict()
    inchikey_category_record = dict()
    mmol_record = dict()
    for k in keys:
        v = record[k]
        if classify_column(k) == COL_REAGENT:
            reagent_record[k] = v
        elif classify_column(k) == COL_REACTION:
            reaction_record[k] = v
        elif classify_column(k) == COL_FEATURE:
            feature_record[k] = v
        elif classify_column(k) == COL_RAW_INCHIKEY_CATEGORY:
            inchikey_category_record[k] = v
        elif classify_column(k) == COL_RAW_MMOL_INCHI_KEY:
            mmol_record[k] = v
    return reagent_record, reaction_record, feature_record, inchikey_category_record, mmol_record


def get_category_x_to_inchikey(inchi_key_type_record: dict):
    # find category for inchi_keys
    category_to_inchikey = dict()
    for k, v in inchi_key_type_record.items():
        if k == "_rxn_organic-inchikey":
            kk = "_raw_organic_0_inchikey"
        else:
            kk = k.lower()
        category = re.sub(r"^_raw_", "", kk)
        category = re.sub(r"_inchikey$", "", category)
        if pd.isna(v):
            continue
        if not isinstance(v, str):
            continue
        category_to_inchikey[category] = v.upper()
    return category_to_inchikey


def get_reagent_table(reagent_record: dict, inventory: list[Material]):
    reagent_dict = dict()
    for i in range(_max_num_reagents):
        rm_dict = dict()
        reagent_volume = f"_raw_reagent_{i}_volume"
        reagent_volume_prepare = f"_raw_reagent_{i}_instructions_2_volume"
        reagent_volume_prepare_unit = f"_raw_reagent_{i}_instructions_2_volume_units"
        try:
            reagent_volume = reagent_record[reagent_volume]
            assert not pd.isna(reagent_volume)
            assert reagent_volume > 1e-6
        except (KeyError, AssertionError) as e:
            continue

        try:
            reagent_volume_prepare = reagent_record[reagent_volume_prepare]
            reagent_volume_prepare_unit = reagent_record[reagent_volume_prepare_unit]
            assert not pd.isna(reagent_volume_prepare) and reagent_volume_prepare > 1e-6 and not pd.isna(
                reagent_volume_prepare_unit)
        except (KeyError, AssertionError) as e:
            reagent_volume_prepare = None
            reagent_volume_prepare_unit = None
            logger.warning(f"reagent preparation volume is missing for reagent index=={i}")

        for j in range(_max_num_rm):
            reagent_chemical_inchikey = f"_raw_reagent_{i}_chemicals_{j}_inchikey"
            reagent_chemical_amount = f"_raw_reagent_{i}_chemicals_{j}_actual_amount"
            reagent_chemical_amount_unit = f"_raw_reagent_{i}_chemicals_{j}_actual_amount_units"
            try:
                reagent_chemical_inchikey = reagent_record[reagent_chemical_inchikey]
                reagent_chemical_amount = reagent_record[reagent_chemical_amount]
                reagent_chemical_amount_unit = reagent_record[reagent_chemical_amount_unit]
                assert not (
                        pd.isna(reagent_chemical_inchikey)
                        or pd.isna(reagent_chemical_amount)
                        or pd.isna(reagent_chemical_amount_unit)
                )
            except (AssertionError, KeyError) as e:
                continue

            rm = ReagentMaterial(
                material=Material.select_from(reagent_chemical_inchikey, inventory),
                amount=reagent_chemical_amount,
                amount_unit=reagent_chemical_amount_unit,
            )
            rm_dict[j] = rm

        if len(rm_dict) == 0:
            continue

        if reagent_volume_prepare_unit == "milliliter":
            assert reagent_volume_prepare is not None
            reagent_volume_prepare *= 1e-3

        reagent = Reagent(
            reagent_material_table=rm_dict,
            volume_added=reagent_volume * 1e-6,  # uL -> L
            volume_prepare=reagent_volume_prepare,
        )
        reagent_dict[i] = reagent
    return reagent_dict


def get_inchikey_to_features(feature_record: dict, category_x_to_inchikey: dict[str, str]) -> dict[
    str, dict[str, float]]:
    data = dict()
    for k, v in feature_record.items():
        chemtype = re.findall(r"^[a-z]*_\d", k[6:])[0]
        feature_name = k.replace(chemtype, "", 1)
        try:
            inchikey = category_x_to_inchikey[chemtype]
        except KeyError:
            continue
        try:
            data[inchikey].update({feature_name: v})
        except KeyError:
            data[inchikey] = {feature_name: v}
    return data


def classify_column(cname: str) -> str:
    if cname == _target_column:
        return COL_TARGET

    if cname == "name":
        return COL_NAME

    if is_column_useless(cname):
        return COL_USELESS

    if cname.startswith("_feat"):
        return COL_FEATURE

    if cname.startswith("_calc"):
        return COL_CALCULATE

    if cname.startswith("_raw_reagent"):
        return COL_REAGENT

    if cname.startswith("_rxn") and not cname.endswith("-inchikey"):
        return COL_REACTION

    if cname.startswith("_raw_mmol") and cname.count("-") == 2:
        return COL_RAW_MMOL_INCHI_KEY

    elif cname.startswith("_raw_molarity") and cname.count("-") == 2:
        return COL_RAW_MOLARITY_INCHI_KEY

    elif (cname.startswith("_raw") and cname.lower().endswith("inchikey")) or (cname == "_rxn_organic-inchikey"):
        # `_rxn_organic-inchikey` seems to be `_raw_organic_0_inchikey`
        return COL_RAW_INCHIKEY_CATEGORY

    elif cname.startswith("_raw"):
        return COL_RAW_OTHER

    raise ValueError(f"unknown column name: {cname}")


def collect_reactions(csv: FilePath, inventory: list[Material], ignore_absent_outcome=False) -> list[Reaction]:
    df = pd.read_csv(csv, low_memory=False)
    records = df.to_dict(orient="records")
    reactions = []
    for r in records:
        try:
            reaction = record_to_reaction(r, inventory=inventory, keys=df.columns,
                                          ignore_absent_outcome=ignore_absent_outcome)
        except ReactionParserError as e:
            logger.critical(str(e))
            logger.critical("the reaction is DROPPED!")
            continue
        reactions.append(reaction)
    return reactions


def is_column_useless(cname: str) -> bool:
    _useless_columns = (
        "_raw_challengeproblem",
        "_raw_genver",
        "_raw_user_generated_experimentname",
        "_raw_batch_count",
        "_raw_datecreated",
        "_raw_jobserial",
        "dataset",
        "_raw_lab",
        "_raw_labwareid",
        "_raw_operator",
        "_raw_modelname",
        "_raw_notes",
        "_raw_modelname",
        "_raw_vialsite",
        "_raw_timecreated_utc",
        "_raw_participantname",
    )
    if cname in _useless_columns:
        return True

    if "nominal_amount" in cname:
        return True

    if cname.endswith("_date"):
        return True

    return False


def record_to_reaction(record: dict, inventory: list[Material], keys: list[str] = None,
                       ignore_absent_outcome=False) -> Reaction:
    name = record['name']
    logger.info(f">>>>> Working on: {name} <<<<<")
    if keys is None:
        keys = list(record.keys())
    reagent_record, reaction_record, feature_record, inchi_key_type_record, mmol_record = get_subrecords(record, keys)

    # this should be bijective
    category_x_to_inchikey = get_category_x_to_inchikey(inchi_key_type_record)
    inchikey_to_category_x = {v: k for k, v in category_x_to_inchikey.items()}
    assert len(category_x_to_inchikey) == len(inchikey_to_category_x)

    reagent_table = get_reagent_table(reagent_record, inventory)

    # get moles from `_raw_mmol`
    mmol_data = dict()
    for k, v in mmol_record.items():
        inchi_key = [word for word in k.split("_") if word.count("-") == 2]
        assert len(inchi_key) == 1
        inchi_key = inchi_key[0].upper()
        mmol = v
        if mmol < 1e-5:
            continue
        if inchi_key in mmol_data:
            mmol_data[inchi_key] += mmol
        else:
            mmol_data[inchi_key] = mmol

    # get reaction time
    try:
        reaction_time = record['_rxn_reactiontime_s']
        assert not pd.isna(reaction_time)
        assert reaction_time > 1e-5
    except (KeyError, AssertionError) as e:
        emsg = f"weird reaction time: {name}"
        raise ReactionParserError(emsg)

    # get reaction temperature
    try:
        reaction_temperature = record['_rxn_temperature_c']
        assert not pd.isna(reaction_temperature)
        assert reaction_temperature > -273.15
    except (KeyError, AssertionError) as e:
        emsg = f"weird reaction temperature: {name}"
        raise ReactionParserError(emsg)

    # get outcome
    outcome = None
    try:
        outcome = record[_target_column]
        assert not pd.isna(outcome)
        outcome = int(outcome)
    except (KeyError, AssertionError, ValueError) as e:
        emsg = f"invalid outcome: {name}, {outcome}"
        if not ignore_absent_outcome:
            raise ReactionParserError(emsg)

    # get expver
    expver = str(record["_raw_expver"])
    try:
        assert not pd.isna(expver)
    except (AssertionError, ValueError) as e:
        emsg = f"invalid expver: {name}, {expver}"
        raise ReactionParserError(emsg)

    reaction = Reaction(
        identifier=name,
        outcome=outcome,
        reaction_time=reaction_time,
        reaction_temperature=reaction_temperature,
        experiment_version=expver,
        reagent_table=reagent_table,
        raw_=record,
        properties={
            "inchikey_to_features": get_inchikey_to_features(feature_record, category_x_to_inchikey),
            "category_x_to_inchikey": category_x_to_inchikey,
            "mmol_data": mmol_data
        }
    )

    # check the obtained reaction
    if reaction.total_volume < 1e-5:
        raise ReactionParserError(f"invalid total volume: {reaction.total_volume}")
    return reaction
