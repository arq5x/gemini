"""
After decompose has split (and before normalize), we still need to adjust the CLN* fields based
on the CLNALLE field.
If we have ALT=A,C,T but CLNALLE=1,3 then there will be only 2 fields in the
rest of the CLN* fields. So we always take the first to be set with ALT=A and
the 2nd to be set with ALT=T and we set C to empty.
"""
import sys
from collections import OrderedDict as dict

for line in sys.stdin:
    if line.startswith('#'):
        print line,
        continue

        print line,
        continue

    fields = line.rstrip('\n').split('\t')
    info = dict((kv[0], (kv[1] if len(kv) > 1 else None)) for kv in (x.split('=') for x in fields[7].split(';')))
    ref, alt = fields[3:5]

    if 'OLD_MULTIALLELIC' in info:
        om = info['OLD_MULTIALLELIC'] # 1:984971:G/A/C
        alleles = om.rsplit(':', 1)[-1].split('/')
        assert ref == alleles[0]
    else:
        if ',' in info['CLNALLE']:
            idxs = map(int, info['CLNALLE'].split(','))

            # make list of alleles like we got from vt decompose so we can use
            # the same logic below.
            alleles = [[ref, alt][i] for i in idxs]
        else:
            print line,
            continue


    # becomes 1-based because alleles[0] is reference.
    alt_i = alleles.index(alt)

    # NOW match this with the order in alleles.
    allele_idxs = map(int, info['CLNALLE'].split(','))
    # special-case the reference allele
    info0 = {}
    for pull_idx, ai in enumerate(allele_idxs):
        if ai == 0:
            for k in (key for key in info if key.startswith('CLN')):
                val = info[k].split(',')
                # in-case we have 0,0 somehow
                if k in info0: info0[k] += ','
                info0[k] = val[pull_idx]

    # there was no CLN record for this alternate
    if not alt_i in allele_idxs:
        for k in (key for key in info if key.startswith('CLN')):
            info[k] = "."
    else:
        # now, allele_idxs will be something like 1, 3.
        # get the correct index for this variant.
        pull_idx = allele_idxs.index(alt_i)
        assert pull_idx >= 0, (allele_idxs, alt_i)
        for k in (key for key in info if key.startswith('CLN')):
            val = info[k].split(',')
            assert len(val) == len(allele_idxs), (info, k, info[k])
            info[k] = val[pull_idx]

        vv = info['CLNALLE'].split(',')
        # adjust the number because now we have only 1 alt.
        info['CLNALLE'] = ",".join(v if v == '0' else '1' for v in vv)

    # add back in the stuff for the reference allele if present.
    for k in info0:
        if info[k] == ".": info[k] = info0[k]
        else: info[k] = info0[k] + "," + info[k]
        

    fields[7] = ";".join(kv[0] + (('=' + kv[1]) if kv[1] is not None else '') for kv in info.items())
    print "\t".join(fields)
