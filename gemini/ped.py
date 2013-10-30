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
        header = possible_header.replace("#", "").split()
        # rename the standard fields to a common name
        header = default_ped_fields + header[len(default_ped_fields):]
        return possible_header.replace("#", "").split()
    else:
        return default_ped_fields

def load_ped_file(ped_file):
    ped_dict = {}
    header = get_ped_fields(ped_file)
    for line in open(ped_file, 'r'):
        if line.startswith("#"):
            continue
        fields = line.split()
        ped_dict[fields[1]] = fields
    return ped_dict
