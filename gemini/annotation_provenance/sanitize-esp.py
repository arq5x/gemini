import sys
from collections import OrderedDict

j = 0
for line in sys.stdin:
    if line[0] == "#":
        if line.startswith(("##INFO=<ID=EA_AC,", "##INFO=<ID=AA_AC,", "##INFO=<ID=TAC,")):

            line = line.replace(",Number=.,Type=String,", ",Number=R,Type=Integer,")
            line = line.replace("AltAlleles,RefAllele", "RefAllele,AltAlleles")
            assert "Number=R" in line
        elif line.startswith("##INFO=<ID=GTS,"):
            line = line.replace(",Number=.", ",Number=G")
            assert "Number=G" in line
        print line,
        continue
    j += 1
    # in body, need to adjust GTS, TAC, AA_AC, EA_AC since ESP stores ref last
    # and Number=R means ref should be first.
    fields = line.split("\t")
    info = OrderedDict((p[0], p[1]) for p in (kv.split("=") for kv in fields[7].split(";")))

    # ignore this for now. need to fix if we use them in the db
    #order = info['GTS'].split(",")
    #A1A1,A1A2,A1A3,A1R,A2A2,A2A3,A2R,A3A3,A3R,RR


    alts = fields[4].split(",")
    for field in ("TAC", "AA_AC", "EA_AC"):
        vals = info[field].split(",")
        assert len(vals) == len(alts) + 1, (fields, j)

        vals = vals[-1:] + vals[:-1]
        assert len(vals) == len(alts) + 1, (fields, j)
        info[field] = ",".join(vals)

    k = len(alts) + 1
    for field in ("GTS", "EA_GTC", "AA_GTC"):

        vals = info[field].split(",")
        if len(vals) != k * (k + 1) / 2:
            # X and Y have incorrect numbers here...
            # need to think about this if we end up using them in db
            # but  the GT fields are currently not used.
            assert fields[0] in "XY"

        #info[field] = ",".join(vals[-1:] + vals[:-1])

    fields[7] = ";".join("%s=%s" % kv for kv in info.items())
    print "\t".join(fields),
