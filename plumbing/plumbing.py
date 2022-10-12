import copy
from collections import defaultdict

from loguru import logger

from escalate_plumber import Material, collect_reactions, WF3Data, json_dump, group_reactions, Reaction, \
    SamplerConvexHull

logger.add(f"{__file__}.log", level="DEBUG")

if __name__ == '__main__':

    ###########################
    # COLLECT REACTIONS
    ###########################
    logger.info("=" * 32)
    logger.info("COLLECT REACTIONS")
    logger.info("=" * 32)
    # build materials inventory
    MaterialInventory = Material.from_csv("../data/ChemicalInventory.csv")
    assert len(MaterialInventory) == len(set(MaterialInventory))

    # collect reactions from csv files, invalid rows will be dropped
    ReactionsIodides = collect_reactions("../data/0064.wf3_iodides.csv", inventory=MaterialInventory,
                                         ignore_absent_outcome=False)
    ReactionsIodidesAlloying = collect_reactions("../data/0064.wf3_iodide_alloying.csv", inventory=MaterialInventory,
                                                 ignore_absent_outcome=False)
    ReactionsAlloys = collect_reactions("../data/0064.wf3_alloys.csv", inventory=MaterialInventory,
                                        ignore_absent_outcome=True)

    set_iodides = set([r.identifier for r in ReactionsIodides])
    set_iodides_alloying = set([r.identifier for r in ReactionsIodidesAlloying])
    set_alloys = set([r.identifier for r in ReactionsAlloys])

    logger.info(f"``iodides` is a proper superset of `iodides alloying`: {set_iodides > set_iodides_alloying}")
    logger.info(f"``iodides` is a proper superset of `alloys`: {set_iodides > set_alloys}")
    logger.info(f"`iodides` intersects `alloys`: {len(set_iodides.intersection(set_alloys))}")
    logger.info(
        f"sizes of `alloys`, `iodides alloying`, `iodides` : {(len(ReactionsAlloys), len(ReactionsIodidesAlloying), len(ReactionsIodides))}")

    # combine reactions, remove duplicates
    Reactions = {r.identifier: r for r in ReactionsIodides + ReactionsIodidesAlloying + ReactionsAlloys if r.is_wf3}
    logger.info("# of wf3 reactions: {}".format(len(Reactions)))

    # collect reactions with a valid outcome score
    """
    Qianxiang:
     in wf3_iodides.csv reactions of 2018-11-02 have 1000 uL formic acid (reagent 5) added (see this cell)
     but the raw files in the Science folder indicate 1000 uL DCM (reagent 6) was added. Which one actually happened?
    Mansoor:
     Based on the ExpDataEntry file, it has to be 1000 uL of DCM (reagent 6), not formic acid (reagent 5).
     I am not sure what caused the error in the CSV file
    """
    ProblematicHeaders = ['2018-11-02', ]
    ReactionsValid = {r.identifier: r for r in Reactions.values() if
                      isinstance(r.outcome, int) and r.experiment_header not in ProblematicHeaders}
    logger.info("# of valid wf3 reactions: {}".format(
        len(ReactionsValid)))  # these are the reactions we are interested

    ###########################
    # COLLECT FEATURES
    ###########################
    logger.info("=" * 32)
    logger.info("COLLECT FEATURES")
    logger.info("=" * 32)

    FeatureDict = {}  # FeatureDict["<inchikey>"]["<feature_name>"] -> v
    all_inchikeys = set()
    for ir, r in ReactionsValid.items():
        r: Reaction
        for ikey in r.inchikeys:
            all_inchikeys.add(ikey)
        inchikey_to_features = r.properties['inchikey_to_features']
        for inchikey, feat_dict in inchikey_to_features.items():
            try:
                FeatureDict[inchikey].update(feat_dict)
            except KeyError:
                FeatureDict[inchikey] = copy.deepcopy(feat_dict)

    inchikey_to_cat = {m.inchikey: m.category for m in MaterialInventory if m.inchikey in all_inchikeys}
    for cat in set(inchikey_to_cat.values()):
        feats = []
        for inchikey, fd in FeatureDict.items():
            if inchikey_to_cat[inchikey] == cat:
                feats.append(tuple(sorted(fd.keys())))
                logger.warning(f"{cat} {inchikey} has # of features: {len(fd)}")
        logger.warning(feats[0])
        assert len(set(feats)) == 1

    ###########################
    # GROUP REACTIONS
    ###########################
    logger.info("=" * 32)
    logger.info("GROUP REACTIONS")
    logger.info("=" * 32)
    # - Reactions with the same `expver` should always be grouped together.
    # - Reactions with the same date (header of the reaction name) should always be grouped together,
    #   we call such a group a `header group`
    # - For each header group, we define its representative reaction that has the largest number of `category_x`,
    #   this representative reaction serves as the template for generating csv
    #   [based on the header, we may be able to find experiment specifications in the `Science` folder]
    # - the fingerprint of a reaction comes from the number of unique organic/inorganic/solvent/acid chemicals
    #   two header groups are combined if their representative reactions share the same fingerprint
    GroupDict = defaultdict(list)
    ConvexHull = None
    for expver, expver_group in group_reactions(list(ReactionsValid.values()),
                                                key_function=lambda x: x.experiment_version).items():
        logger.info(f"expver == {expver}: # of reactions == {len(expver_group)}")
        header_to_group = group_reactions(expver_group, key_function=lambda x: x.experiment_header)
        for header, header_group in header_to_group.items():
            logger.info(f"\theader == {header}: # of reactions == {len(header_group)}")
            # a header group is represented by the reaction with most `category_x`
            representative_reaction = sorted(
                header_group, key=lambda x: len(x.properties["category_x_to_inchikey"].keys()), reverse=True
            )[0]
            representative_reaction: Reaction
            try:
                assert all(
                    set(
                        r.properties["category_x_to_inchikey"].keys()).issubset(
                        set(representative_reaction.properties["category_x_to_inchikey"].keys())
                    ) for r in header_group
                )
            except AssertionError:
                for r in header_group:
                    print(r.properties["category_x_to_inchikey"].keys(), )
                raise ValueError()

            wf3data = WF3Data.from_reaction(representative_reaction)

            reagent_set = set()
            for r in header_group:
                reagent_set_except_antisolvent = set(
                    [
                        reagent for reagent in r.reagent_set if not sorted(
                        reagent.molarity_table.keys()
                    )[0] == wf3data.antisolvent_identity
                    ]
                )
                reagent_set = reagent_set.union(reagent_set_except_antisolvent)
            sch = SamplerConvexHull.from_reagent_set(reagent_set)

            if ConvexHull is None:
                ConvexHull = sch
            else:
                ConvexHull = SamplerConvexHull.combine(ConvexHull, sch)

            group_key = "expver-{}_{}".format(expver, wf3data.fingerprint)
            GroupDict[group_key] += [r.identifier for r in header_group]

    assert ConvexHull is not None
    ###########################
    # CONVERT REACTIONS TO DATA
    ###########################
    logger.info("=" * 32)
    logger.info("CONVERT REACTIONS TO DATA")
    logger.info("=" * 32)
    WF3Entries = {i: WF3Data.from_reaction(r) for i, r in ReactionsValid.items()}

    PlumbingData = {
        "ReactionsValid": ReactionsValid,
        "GroupDict": GroupDict,
        "FeatureDict": FeatureDict,
        "WF3Entries": WF3Entries,
        "MaterialInventory": MaterialInventory,
        "ConvexHull": ConvexHull,
    }
    json_dump(PlumbingData, "PlumbingData.json.gz", gz=True)
