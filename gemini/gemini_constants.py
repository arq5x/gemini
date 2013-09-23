#!/usr/bin/env python

VARIANTS_KEY = "variant_id"
BUFFER_SIZE  = 10000


# genotype encoding.
# 0 / 00000000 hom ref
# 1 / 00000001 het
# 2 / 00000010 unknown
# 3 / 00000011 hom alt
HOM_REF = 0
HET = 1
UNKNOWN = 2
HOM_ALT = 3

MISSING = None
UNAFFECTED = 1
AFFECTED = 2
