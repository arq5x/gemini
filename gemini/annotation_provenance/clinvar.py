from __future__ import print_function
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
        print(line, end="")
        continue

        print(line, end="")
        continue

    fields = line.rstrip('\n').split('\t')
    info = dict((kv[0], (kv[1] if len(kv) > 1 else None)) for kv in (x.split('=') for x in fields[7].split(';')))
    ref, alt = fields[3:5]
    assert not "," in alt, "use a new clinvar version"
    assert not 'CLNALLE' in info, "use a new clinvar version"
    assert not 'OLD_MULTIALLELIC' in info, "use a new clinvar version"
    print(line, end="")
