#!/usr/bin/env python

import sys

default_ped_fields = ["family_id", "name", "paternal_id", "maternal_id",
                      "sex", "phenotype"]

def get_ped_fields(ped_file):
    if not ped_file:
        return default_ped_fields

    with open(ped_file) as in_handle:
        possible_header = in_handle.readline()

    if possible_header.startswith("#"):
        return possible_header.replace("#", "").split()
    else:
        return default_ped_fields
