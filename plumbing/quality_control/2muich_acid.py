from escalate_plumber import json_load, WF3Data, Reaction

PlumbingData = json_load("../PlumbingData.json.gz")
ReactionsValid = PlumbingData["ReactionsValid"]
WF3Entries = PlumbingData["WF3Entries"]

if __name__ == '__main__':

    f = open('2much_acid.csv', 'w')
    print('reaction identifier,formic acid molarity in alpha vial (M)', file=f)
    for i in ReactionsValid:
        r = ReactionsValid[i]
        wf3 = WF3Entries[i]
        r: Reaction
        wf3: WF3Data

        if any(m > 10 for m in wf3.acid.values()):
            print(r.identifier + "," + "{:.4f}".format(sorted(wf3.acid.values())[0]), file=f)
    f.close()
