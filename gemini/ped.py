"""Parse input PED files for loading into GEMINI sample table.

Tries to parse as tab-delimited if files contain tabs, otherwise splits
on whitespace.
"""

default_ped_fields = ["family_id", "name", "paternal_id", "maternal_id",
                      "sex", "phenotype"]
missing_member = set(["None", None, "0", "-9"])

def get_ped_fields(ped_file):
    if not ped_file:
        return default_ped_fields

    with open(ped_file) as in_handle:
        possible_header = in_handle.readline()

    if possible_header.startswith("#"):
        if possible_header.count("\t") > 1:
            header = possible_header.replace("#", "").rstrip().split("\t")
        else:
            header = possible_header.replace("#", "").split()
        # rename the standard fields to a common name
        header = default_ped_fields + header[len(default_ped_fields):]
        return header
    else:
        return default_ped_fields

def load_ped_file(ped_file):
    ped_dict = {}
    for line in open(ped_file, 'r'):
        if line.startswith("#") or len(line.strip()) == 0:
            continue
        parts = line.rstrip().split("\t") if line.count("\t") > 1 else line.split()
        fields = [x.strip() for x in _fix_ped_family_fields(parts)]
        fields = _fix_ped_affected(fields)
        fields = _fix_ped_sex(fields)
        ped_dict[fields[1]] = fields
    return ped_dict

def _fix_ped_sex(fields):
    l = {'male': '1', 'female': '2'}
    fields[4] = l.get(fields[4].lower(), fields[4])
    return fields

def _fix_ped_affected(fields):
    l = {'affected': '2', 'unaffected': '1'}
    fields[5] = l.get(fields[5].lower(), fields[5])
    return fields

def _fix_ped_family_fields(fields):
    """
    translate commonly used values for missing family members to the canonical
    0 == missing style used in the PED spec
    """
    family_fields = [0, 2, 3]
    for field in family_fields:
        if fields[field] in missing_member:
            fields[field] = "0"
    return fields
