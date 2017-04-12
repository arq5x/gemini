import sys

for line in sys.stdin:
    toks = line.rstrip().split("\t")
    toks[3] = toks[3].replace("?", "")
    if toks[3].startswith("DNA"):
        toks[3] = "DNAtr"
    elif toks[3] == "Low_complexity_Low_complexity":
        toks[3] = "LCR"
    elif toks[3].startswith("LTR"):
        toks[3] = "LTR"
    elif "satellite" in toks[3].lower():
        toks[3] = "satellite"
    elif  toks[3].startswith("Simple_repeat"):
        toks[3] = "simple"
    elif toks[3].split("_")[0] == toks[3].split("_")[1]:
        toks[3] = toks[3].split("_")[0].lower()
    elif toks[3].startswith("LINE_") and toks[3] not in ("LINE_L1", "LINE_L2"):
        toks[3] = "LINE"
    elif toks[3].startswith("SINE_") and toks[3] != "SINE_Alu":
        toks[3] = "SINE"
    if toks[3] == "sine":
        toks[3] = "SINE"
    if toks[3] == "other":
        toks[3] = "unknown"
    print("\t".join(toks))
