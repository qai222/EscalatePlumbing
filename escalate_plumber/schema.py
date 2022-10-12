from __future__ import annotations

import itertools
import re
import typing

import numpy as np
import pandas as pd
from loguru import logger
from monty.json import MSONable
from pandas._typing import FilePath

EscalateCategories = ['organic', 'inorganic', 'solvent', 'acid']


class Material(MSONable):
    def __init__(self, inchikey: str, mw: float, category: str, name: str, density: float, inchi: str,
                 properties: dict = None):
        """
        describe a chemical

        :param inchikey:
        :param mw: in g/mol
        :param category:
        :param name:
        :param density: in g/mL
        :param properties:
        """
        self.inchi = inchi
        self.properties = properties
        self.density = density
        self.name = name
        self.category = category
        self.mw = mw
        self.inchikey = inchikey

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.inchikey)

    def __lt__(self, other: Material):
        return self.inchikey.__lt__(other.inchikey)

    def __gt__(self, other: Material):
        return self.inchikey.__gt__(other.inchikey)

    @property
    def pure_molarity(self):
        return 1e3 * self.density / self.mw

    @staticmethod
    def select_from(inchikey: str, mat_list: list[Material]):
        for m in mat_list:
            if inchikey == m.inchikey:
                return m
        raise ValueError(f"Material not found: {inchikey}")

    @staticmethod
    def from_csv(csv: FilePath) -> list[Material]:
        dataframe = pd.read_csv(csv, low_memory=False)
        mats = []
        for r in dataframe.to_dict(orient="records"):
            inchikey = r["InChI Key (ID)"]
            inchi = r["InChI="]
            name = r["Chemical Name"]
            mw = r["Molecular Weight (g/mol)"]
            density = r["Density            (g/mL)"]
            category_string = r["Chemical Category"]
            if any(pd.isna(v) for v in [inchikey, name, mw, density, category_string]):
                continue
            mw = float(str(mw).replace(",", ""))
            density = float(density)
            category_string = category_string.strip().lower()
            if "," in category_string:
                cts = category_string.split(",")
                if "acid" in cts:
                    category = 'acid'
                else:
                    raise ValueError(f"more than one category: {category_string}")
            else:
                category = category_string
            mat = Material(
                inchikey=inchikey, mw=mw, category=category, name=name, density=density, inchi=inchi
            )
            mats.append(mat)
        return mats


class ReagentMaterial(MSONable):

    def __init__(self, material: Material, amount: float, amount_unit: str):
        self.material = material
        self.amount = amount
        self.amount_unit = amount_unit
        assert self.amount_unit.lower() in ("milliliter", "gram"), f"unconventional amount unit: {self.amount_unit}"

    @property
    def volume(self) -> float:  # in Liter
        if self.amount_unit.lower() == "gram":
            v = 1e-3 * self.amount / self.material.density
        elif self.amount_unit.lower() == "milliliter":
            v = 1e-3 * self.amount
        else:
            raise ValueError(f"unknown amount unit: {self.amount_unit}")
        return v

    @property
    def mol(self) -> float:
        if self.amount_unit.lower() == "gram":
            gram = self.amount
        elif self.amount_unit.lower() == "milliliter":
            gram = self.material.density * self.amount
        else:
            raise ValueError(f"unknown amount unit: {self.amount_unit}")
        return gram / self.material.mw

    def __repr__(self):
        return f"{self.__class__.__name__} {self.material}: {self.amount} {self.amount_unit}"


class Reagent(MSONable):

    def __init__(self, reagent_material_table: dict[int, ReagentMaterial],
                 volume_added: float, volume_sum: float = None,
                 molarity_table: dict[str, float] = None):
        """
        describing a reagent added directly to the reaction, can be a mixture

        Joshua: measured volume in `instruction` maybe grossly inaccurate, better use volume sum

        :param reagent_material_table:
        :param volume_added: in Liter
        :param molarity_table:
        """
        self._volume_sum = volume_sum
        self.volume_added = volume_added
        self.reagent_material_table = reagent_material_table
        self._molarity_table = molarity_table

        assert len(self.reagent_material_table) == len(
            set([c.material for c in self.reagent_material_table.values()])
        ), "duplicate materials in the reagent!"

    @property
    def volume_sum(self):
        """ Joshua: volume is always additive, so it can be calculated from reagent materials """
        vol_sum = 0
        for c in self.reagent_material_table.values():
            vol_sum += c.volume
        return vol_sum

    @property
    def molarity_table(self) -> dict[str, float]:
        d = dict()
        for ic, c in self.reagent_material_table.items():
            d[c.material.inchikey] = c.mol / self.volume_sum
        return d

    def __repr__(self):
        return "\n".join(
            [f"{self.__class__.__name__}", f"\t\tadded volume {self.volume_added} Liter:",
             f"\t\tstock solution volume from volume sum: {self.volume_sum} Liter"
             ] +
            [f"\t\t i=={i}: {c}" for i, c in self.reagent_material_table.items()]
        )

    def __eq__(self, other: Reagent):
        if not isinstance(other, Reagent):
            return False
        if len(self.molarity_table) != len(other.molarity_table):
            return False
        if set(self.molarity_table.keys()) != set(other.molarity_table.keys()):
            return False
        for k in self.molarity_table.keys():
            if abs(self.molarity_table[k] - other.molarity_table[k]) > 1e-5:
                return False
        return True

    def __hash__(self):
        return hash(tuple(sorted([(k, round(v, 5)) for k, v in self.molarity_table.items()])))


class Reaction(MSONable):

    def __init__(
            self,
            identifier: str,
            outcome: int,

            reaction_time: float,
            reaction_temperature: float,
            experiment_version: str,

            reagent_table: dict[int, Reagent],

            raw_: dict = None,
            properties: dict = None,
    ):
        self.reagent_table = reagent_table
        self.experiment_version = experiment_version
        self.reaction_time = reaction_time
        self.reaction_temperature = reaction_temperature
        self.outcome = outcome
        self.identifier = identifier

        if raw_ is None:
            raw_ = dict()
        if properties is None:
            properties = dict()
        self.properties = properties
        self.raw_ = raw_

    @property
    def reagent_set(self):
        return set([r for r in self.reagent_table.values()])

    @property
    def experiment_header(self):
        date_regex_minute = r"^\d{4}\-\d{2}-\d{2}T\d{2}_\d{2}"
        date_regex_day = r"^\d{4}\-\d{2}-\d{2}"
        try:
            date = re.findall(pattern=date_regex_minute, string=self.identifier)
            assert len(date) > 0
        except AssertionError:
            date = re.findall(pattern=date_regex_day, string=self.identifier)
        assert len(date) == 1
        return date[0]

    @property
    def total_volume(self):
        return sum(r.volume_added for r in self.reagent_table.values())  # in Liter

    def __repr__(self):
        return "\n".join(
            [f"{self.__class__.__name__}: {self.identifier}", ] +
            [f"\ti=={i}: {c}" for i, c in self.reagent_table.items()]
        )

    @property
    def inchikey_to_material(self) -> dict[str, Material]:
        d = dict()
        for r in self.reagent_table.values():
            for rm in r.reagent_material_table.values():
                d[rm.material.inchikey] = rm.material
        return d

    @property
    def inchikeys(self):
        return sorted(self.inchikey_to_material.keys())

    @property
    def is_wf3(self) -> bool:
        return self.experiment_version.startswith("3")


class WF3Data(MSONable):

    def __init__(
            self,
            identifier: str,
            fingerprint: str,
            outcome: int,
            # molarities,
            organic: dict[str, float],
            inorganic: dict[str, float],
            solvent: dict[str, float],
            acid: dict[str, float],
            # volume of the reactor,
            alpha_vial_volume: float,
            # volume of the antisolvent,
            beta_vial_volume: float,
            # reaction parameters,
            reaction_time: float,
            reaction_temperature: float,
            antisolvent_identity: str,
    ):

        self.identifier = identifier
        self.fingerprint = fingerprint
        self.outcome = outcome
        self.organic = organic
        self.inorganic = inorganic
        self.solvent = solvent
        self.acid = acid
        self.alpha_vial_volume = alpha_vial_volume
        self.beta_vial_volume = beta_vial_volume
        self.reaction_time = reaction_time
        self.reaction_temperature = reaction_temperature
        self.antisolvent_identity = antisolvent_identity

    @classmethod
    def from_reaction(cls, reaction: Reaction):
        logger.info(">" * 12 + f"CONVERT TO `WF3Data` FROM REACTION: {reaction.identifier}")
        assert reaction.is_wf3

        # a wf3 reaction must have one and only one antisolvent
        # there should be one reagent contains only one chemical,
        # and this reagent was added in a significant amount (>500 uL)
        possible_antisolvents = set()
        for ireagent, reagent in reaction.reagent_table.items():
            if len(reagent.reagent_material_table) == 1 and reagent.volume_added > 500 * 1e-6:
                possible_antisolvents.add(ireagent)
        assert len(possible_antisolvents) == 1
        ireagent_antisolvent = list(possible_antisolvents)[0]
        reagent_antisolvent = reaction.reagent_table[ireagent_antisolvent]
        ireagent_material_antisolvent = list(reagent_antisolvent.reagent_material_table.keys())[0]
        reagent_material_antisolvent = reagent_antisolvent.reagent_material_table[ireagent_material_antisolvent]

        # vial volumes
        volume_alpha_vial = 0  # in liter
        volume_beta_vial = 0  # in liter
        for ireagent, reagent in reaction.reagent_table.items():
            if ireagent == ireagent_antisolvent:
                volume_beta_vial += reagent.volume_added
            else:
                volume_alpha_vial += reagent.volume_added

        # molarities
        # there are at least 4 possible ways to get molarities/mole amounts added to the system:
        #
        # 1. use _raw_molarity_<inchikey> directly, for workflow 3 this is problematic: it was calculated using alpha + beta volume
        #
        # 2. use _raw_molarity_<inchikey> but convert it back to mol using total volume:
        # volume of alpha vial (sum of all _raw_reagent_x_volume that is not DCM) +
        # volume of beta vial (_raw_reagent_x_volume of DCM),
        # then calculate molarities using the volumes of individual vials
        #
        # 3. use _raw_mmol_<inchikey> and the volumes in 2.
        # this has been confirmed to be equivalent to 2.
        #
        # 4. start from _raw_reagent_i_chemicals_j_actual_amount,
        # use molecular weight (and density if liquid) to calculate the mole amount of chemical_j in reagent_i,
        # then use the volume of prepared reagent_i to get the molarity of chemical_j in reagent_i.
        # This molarity times the volume of reagent_i actually added to the system is the mol amount of chemical_j.

        use_mmol_data = "mmol_data" in reaction.properties
        use_reagent_data = True
        assert use_mmol_data

        category_molarity_table_mmol_data = {cat: {} for cat in EscalateCategories}
        category_molarity_table_reagent_data = {cat: {} for cat in EscalateCategories}
        if use_mmol_data:
            logger.warning("using `_raw_mmol_<inchikey>` to calculate molarity")
            mmol_data = reaction.properties['mmol_data']
            for inchikey, mat in reaction.inchikey_to_material.items():
                mmol = mmol_data[inchikey]
                if inchikey == reagent_material_antisolvent.material.inchikey:
                    vol = volume_beta_vial
                else:
                    vol = volume_alpha_vial
                molarity = mmol * 1e-3 / vol
                category_molarity_table_mmol_data[mat.category][inchikey] = molarity

                logger.info(f"inchikey {inchikey} ({mat.name}) has mol amount: {mmol * 1e-3} mol")
                logger.info(f"the volume is: {vol} Liter")
                logger.info(f"the molarity is: {molarity}")

        if use_reagent_data:
            logger.warning("using individual reagent molarities to calculate molarity")
            for ireagent, reagent in reaction.reagent_table.items():
                assert reagent.molarity_table is not None
                for irm, rm in reagent.reagent_material_table.items():
                    moles = reagent.molarity_table[rm.material.inchikey] * reagent.volume_added
                    if ireagent == ireagent_antisolvent:
                        vol = volume_beta_vial
                    else:
                        vol = volume_alpha_vial
                    try:
                        category_molarity_table_reagent_data[rm.material.category][rm.material.inchikey] += moles / vol
                    except KeyError:
                        category_molarity_table_reagent_data[rm.material.category][rm.material.inchikey] = moles / vol

                    logger.info(
                        f"inchikey {rm.material.inchikey} ({rm.material.name}) in ireagent=={ireagent} has mol amount: {moles} mol")
                    logger.info(f"the volume is: {vol} Liter")
                    logger.info(f"the molarity is: {moles / vol}")

        if use_reagent_data and use_mmol_data:
            logger.info("both mmol_data and reagent_data present, check if they are similar")
            assert set(category_molarity_table_reagent_data.keys()) == set(category_molarity_table_mmol_data.keys())
            for inchikey, mat in reaction.inchikey_to_material.items():
                molarity_from_mmol = category_molarity_table_mmol_data[mat.category][inchikey]
                molarity_from_reagent = category_molarity_table_reagent_data[mat.category][inchikey]

                logger.info(f"molarity of {inchikey} ({mat.name}) from reagent_data: {molarity_from_reagent} M")
                logger.info(f"molarity of {inchikey} ({mat.name}) from mmol_data: {molarity_from_mmol} M")
                adiff = abs(molarity_from_reagent - molarity_from_mmol)
                rdiff = adiff / molarity_from_reagent
                if rdiff < 1e-3:
                    logger.info(f"the difference is {adiff} M {rdiff}%")
                else:
                    logger.critical(f"the difference is {adiff} M {rdiff}% ??????")
                # this difference may come from the fact that `mmol_data` was (maybe) calculated using nominal amount

        del category_molarity_table_mmol_data[
            reagent_material_antisolvent.material.category
        ][
            reagent_material_antisolvent.material.inchikey
        ]

        fingerprint = "%".join([str(len(category_molarity_table_mmol_data[k])) for k in EscalateCategories])

        return cls(
            identifier=reaction.identifier,
            fingerprint=fingerprint,
            outcome=reaction.outcome,
            alpha_vial_volume=volume_alpha_vial,
            beta_vial_volume=volume_beta_vial,
            reaction_time=reaction.reaction_time,
            reaction_temperature=reaction.reaction_temperature,
            antisolvent_identity=reagent_material_antisolvent.material.inchikey,
            **category_molarity_table_mmol_data,
        )


class SamplerConvexHull(MSONable):
    def __init__(self, space: typing.Tuple[Material, ...], df: pd.DataFrame):
        self._space = space
        self._df = df

    @staticmethod
    def combine(sch1: SamplerConvexHull, sch2: SamplerConvexHull) -> SamplerConvexHull:
        new_space = []
        for mat in sch1.space:
            new_space.append(mat)
        for mat in sch2.space:
            new_space.append(mat)
        new_space = tuple(sorted(set(new_space)))
        records = sch1.df.to_dict(orient='records') + sch2.df.to_dict(orient='records')
        new_df = pd.DataFrame.from_records(records, columns=[m.inchikey for m in new_space])
        new_df.fillna(value=0.0, inplace=True)
        new_df.drop_duplicates(inplace=True, ignore_index=True)
        return SamplerConvexHull(new_space, new_df)

    @property
    def space(self):
        assert len(self._space) > 0
        return tuple(sorted(self._space))

    @property
    def df(self):
        cols = [m.inchikey for m in self.space]
        assert cols == [c for c in self._df.columns]
        return self._df.sort_values(by=cols, axis=0)

    def __eq__(self, other: SamplerConvexHull):
        if len(self.space) != len(other.space):
            return False
        if any(self.space[i] != other.space[i] for i in range(len(self.space))):
            return False
        if self.df.shape != other.df.shape:
            return False
        if not np.allclose(self.df.values, other.df.values):
            return False
        return True

    def __repr__(self):
        s = self.__class__.__name__
        s += f"\n# of axes: {len(self.space)}"
        s += f"\naxes:"
        for mat in self.space:
            s += f"\n\t{mat.name} {mat.inchikey}"
        s += "\n"
        s += self.df.__repr__()
        return s

    @classmethod
    def from_reagent_set(cls, reagent_set: set[Reagent]):
        concentration_space = set()
        for r in reagent_set:
            for rmat in r.reagent_material_table.values():
                concentration_space.add(rmat.material)
        concentration_space = tuple(sorted(concentration_space))
        points = []
        for r in reagent_set:
            v = [0., ] * len(concentration_space)
            for iax, ax in enumerate(concentration_space):
                try:
                    v[iax] = r.molarity_table[ax.inchikey]
                except KeyError:
                    v[iax] = 0.
            assert len(v) > 0 and any(vv > 1e-7 for vv in v)
            points.append(v)
        df = pd.DataFrame(points, columns=[m.inchikey for m in concentration_space])
        return cls(concentration_space, df)


def group_reactions(reactions: typing.Union[list[Reaction], list[WF3Data]], key_function: typing.Callable):
    groups = []
    unique_keys = []
    data = sorted(reactions, key=key_function)
    for k, g in itertools.groupby(data, key_function):
        groups.append(list(g))
        unique_keys.append(k)
    return dict(zip(unique_keys, groups))
