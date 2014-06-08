import sys
import snpEff
import vep


def interpret_impact(args, var, effect_fields):
    """
    Interpret the report from SnpEff or VEP to determine the impact of the variant.

    SnpEff examples:
    0    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008),
    1    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000534990|exon_1_69037_69829)

    VEP examples
    CSQ: Consequence|Codons|Amino_acids|Gene|hgnc|Feature|EXON|polyphen|sift
    non_synonymous_codon|gaT/gaG|D/E|ENSG00000116254|CHD5|ENST00000378006|18/25|benign(0.011)|tolerated(0.3)
    nc_transcript_variant|||ENSG00000116254|CHD5|ENST00000491020|5/6||
    """
    impact_all = []
        # holds a list of all the transcript impacts for this variant
    effect_strings_str = ""
    effect_strings = []
    counter = 0  # counter for anno_id
    if args.anno_type == "snpEff":
        try:
            effect_strings_str = var.INFO["EFF"]
            effect_strings = effect_strings_str.split(",")
        except KeyError:
            if "SNPEFF_EFFECT" in var.INFO:
                impact_all.append(snpEff.gatk_effect_details(var.INFO))
            else:
                sys.stderr.write("WARNING: The input VCF has no snpEFF annotations. "
                                 "Variant impact will be set to unknown\n")

        for effect_string in effect_strings:
            counter += 1
            eff_pieces = snpEff.eff_search.findall(effect_string)
            for piece in eff_pieces:
                impact_string = piece[0]
                    # the predicted inpact, which is outside the ()
                impact_detail = piece[1]
                    # all the other information, which is inside the ()
                try:
                    impact_info = snpEff.effect_map[impact_string]
                    impact_details = snpEff.EffectDetails(impact_string,
                                                      impact_info.priority,
                                                      impact_detail,
                                                      counter,
                                                      args.maj_version)
                except KeyError:
                    impact_details = snpEff.EffectDetails(impact_string,
                                                    None,
                                                    impact_detail,
                                                    counter,
                                                    args.maj_version)
                impact_all.append(impact_details)

    elif args.anno_type == "VEP":
        try:
            effect_strings_str = var.INFO["CSQ"]
            effect_strings = effect_strings_str.split(",")
        except KeyError:
            sys.stderr.write("WARNING: The input VCF has no VEP annotations. \
                             Variant impact will be set to unknown\n")

        for effect_string in effect_strings:
             # nc_transcript_variant&intron_variant|||ENSG00000243485|MIR1302-11|ENST00000
            each_string = effect_string.split("|")
            if "&" in each_string[0]:
                impact_strings = each_string[0].split("&")
                # impact_strings will be [nc_transcript_variant,
                # intron_variant]
                for impact_string in impact_strings:
                    counter += 1
                    try:
                        impact_info = vep.effect_map[impact_string]
                        impact_details = vep.EffectDetails(
                            impact_string, impact_info.priority, effect_string, counter,
                            effect_fields)
                    except KeyError:
                        impact_details = vep.EffectDetails(
                            impact_string, None, effect_string, counter, effect_fields)
                    impact_all.append(impact_details)
            # we expect VEP to produce a valid impact label for each_string[0]
            elif "&" not in each_string[0]:
                counter += 1
                impact_string = each_string[0]
                impact_info = vep.effect_map.get(impact_string)
                try:
                    impact_details = vep.EffectDetails(
                        impact_string, impact_info.priority, effect_string, counter, effect_fields)
                except AttributeError:
                    impact_details = vep.EffectDetails(
                        impact_string, None, effect_string, counter, effect_fields)
                impact_all.append(impact_details)
    else:
        # should not get here, as the valid -t options should be handled
        # in main()
        sys.exit("ERROR: Unsupported variant annotation type.\n")

    return impact_all
