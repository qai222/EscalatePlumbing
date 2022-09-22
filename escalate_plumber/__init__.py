"""
useful docs on gdrive

chemical inventory
https://docs.google.com/spreadsheets/d/1JgRKUH_ie87KAXsC-fRYEw_5SepjOgVt7njjQBETxEg

All runs of morph phase mapping papers
https://docs.google.com/document/d/1lJ-xvI2OVibl-BagMEYoYMD9KHaP5PbP/

Science folder
https://drive.google.com/drive/folders/1BMF9imJakX3MqRKb4ou3F0i5W-XAV2-1
"""
import gzip
import json

from monty.json import MontyEncoder, MontyDecoder
from pandas._typing import FilePath

from .schema import Material, Reaction, Reagent, ReagentMaterial, group_reactions, WF3Data, EscalateCategories
from .parser import collect_reactions, classify_column, is_column_useless, record_to_reaction


def json_load(fn: FilePath):
    if fn.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open
    with fopen(fn, "rt") as f:
        return json.load(f, cls=MontyDecoder)


def json_dump(o, fn: FilePath, gz=True):
    if gz:
        fopen = gzip.open
    else:
        fopen = open
    with fopen(fn, 'wt') as f:
        json.dump(o, f, cls=MontyEncoder)
