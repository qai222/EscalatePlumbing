import pprint

import numpy as np
import pandas as pd
from loguru import logger
from tqdm import tqdm

from escalate_plumber import json_load, Reaction, Reagent

PlumbingData = json_load("../PlumbingData.json.gz")
ReactionsValid = PlumbingData["ReactionsValid"]
ConvexHull = PlumbingData["ConvexHull"]

DfHull = pd.read_csv("../csv/convexhull.csv")


def is_reagent_antisolvent(reagent: Reagent):
    """ DCM always functions as the antisolvent """
    return len(reagent.molarity_table) == 1 and \
           sorted(reagent.molarity_table.keys())[0] == "YMWUJEATGCHHMB-UHFFFAOYSA-N"


def drop_zeros(d: dict[str, float], eps=1e-9):
    """ drop dummy entries for molarity table """
    return {k: v for k, v in d.items() if abs(v) > eps}


def get_stock_solution_row_index(molarity_table: dict[str, float], df_hull: pd.DataFrame, eps=1e-7):
    """ given a molarity table, find its row index in `convexhull.csv` """
    mt = drop_zeros(molarity_table)
    found = None
    for i, row in df_hull.iterrows():
        row_no_zeros = drop_zeros(row)
        if set(mt.keys()) != set(row_no_zeros.keys()):
            continue
        if any(mt[k] - row[k] > eps for k in mt):
            continue
        else:
            found = i
            break
    try:
        assert found is not None
    except AssertionError:
        emsg = f"cannot find index for this stock solution: {pprint.pformat(mt)}"
        logger.critical(emsg)
        raise RuntimeError("impossible...")
    return found


def get_reaction_molarity_table(r: Reaction) -> dict[str, float]:
    """ get the molarity table for chemicals in alpha vial """
    mt = dict()
    vol_alpha = 0.0

    for rr in r.reagent_table.values():
        if is_reagent_antisolvent(rr):
            continue
        vol_alpha += rr.volume_added

    for rr in r.reagent_table.values():
        if is_reagent_antisolvent(rr):
            continue
        for i, m in rr.molarity_table.items():
            try:
                mt[i] += m * rr.volume_added / vol_alpha
            except KeyError:
                mt[i] = m * rr.volume_added / vol_alpha
    return mt


def get_checking_params(r: Reaction):
    logger.info(f"checking reaction: {r.identifier}")

    reaction_molarity_table = get_reaction_molarity_table(r)
    inchikeys = sorted(reaction_molarity_table.keys())
    logger.info(f"chemicals: {inchikeys}")

    target_molarity = [reaction_molarity_table[i] for i in inchikeys]

    reagents = [rr for rr in r.reagent_table.values() if not is_reagent_antisolvent(rr)]

    ss_matrix = np.zeros((len(reagents), len(inchikeys)))
    ss_vols = []
    ss_indices = []
    for irr, rr in enumerate(reagents):
        rr: Reagent
        ss_indices.append(get_stock_solution_row_index(rr.molarity_table, DfHull))
        ss_vols.append(rr.volume_added)
        for ii, inchik in enumerate(inchikeys):
            try:
                ss_matrix[irr][ii] = rr.molarity_table[inchik]
            except KeyError:
                continue
    logger.info(f"stock solution molarity matrix, S:\n{ss_matrix}")
    logger.info(f"added volumes, V:\n{ss_vols}")
    logger.info(f"row indices in `convexhull.csv`:\n{ss_indices}")
    logger.info(f"final alpha vial molarities: {target_molarity}")
    return dict(
        target_molarity=target_molarity,
        ss_matrix=ss_matrix,
        ss_vols=ss_vols,
        ss_indices=ss_indices,
        inchikeys=inchikeys,
    )


def forward_check(params: dict):
    """ VS = target * sum(S) """
    ss_vols = params['ss_vols']
    ss_matrix = params['ss_matrix']
    target = params['target_molarity']
    vs = ss_vols @ ss_matrix / sum(ss_vols)
    logger.info(f"VS/alpha_volume: {vs}")
    if np.allclose(vs, target):
        logger.info("FORWARD CHECK: PASSED")
        return True
    else:
        logger.critical("FORWARD CHECK: FAILED")
        return False


def backward_check(target_molarities: dict[str, float], alpha_volume: float, df_hull: pd.DataFrame):
    """ solve volume values as an LP using gurobi """

    import gurobipy as gp
    from gurobipy import GRB

    inchikeys = [str(c) for c in df_hull.columns]
    assert set(inchikeys).issuperset(set(target_molarities.keys()))
    target_molarities_padded = np.zeros(len(df_hull.columns))
    for i, inchikey in enumerate(inchikeys):
        try:
            target_molarities_padded[i] = target_molarities[inchikey]
        except KeyError:
            continue

    n_ss, n_chem = df_hull.shape

    # init gurobi model, suppress output
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)
        env.setParam('LogToConsole', 0)
        env.start()
        with gp.Model(env=env) as m:

            v = []
            for i in range(n_ss):
                v_i = m.addVar(lb=0., name=str(i))  # add var
                v.append(v_i)
                m.addConstr(v_i >= 0., name=str(i))

            m.addConstr(sum(v) == alpha_volume, name='volume_sum')
            for j in range(n_chem):
                c = 0.
                for i in range(n_ss):
                    c += v[i] * df_hull.values[i][j]
                m.addConstr(c == alpha_volume * target_molarities_padded[j], name=inchikeys[j])

            # objective function
            objective = 1

            m.setObjective(objective, GRB.MINIMIZE)
            m.optimize()

            try:
                sol = [v.x for v in m.getVars()]
                logger.info("BACKWORD PASSED")
                return True
            except AttributeError:
                logger.critical("BACKWORD FAILED")
                return False


if __name__ == '__main__':
    logger.remove()
    logger.add(f"{__file__}.log", level="DEBUG")

    for r in tqdm(ReactionsValid.values()):
        params = get_checking_params(r)
        t_mol = dict(zip(params['inchikeys'], params['target_molarity']))
        a_vol = sum(params['ss_vols'])
        assert forward_check(params)
        assert backward_check(t_mol, a_vol, DfHull)
