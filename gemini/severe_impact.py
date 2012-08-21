import sys
import snpEff, vep

def interpret_severe_impact(args, var):
    """
    Interpret the report from SnpEff or VEP to determine the impact of the variant.
    
    SnpEff examples:
    0    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008),
    1    NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T/A|OR4F5|protein_coding|CODING|ENST00000534990|exon_1_69037_69829)
    
    VEP examples
    CSQ: Consequence|Codons|Amino_acids|Gene|hgnc|Feature|EXON|polyphen|sift|condel
    non_synonymous_codon|gaT/gaG|D/E|ENSG00000116254|CHD5|ENST00000378006|18/25|benign(0.011)|tolerated(0.3)|neutral(0.029)
    nc_transcript_variant|||ENSG00000116254|CHD5|ENST00000491020|5/6|||
    """
    
    effect_strings_str = ""
    effect_strings = []
    impact_details = None
    max_severity = 9 # initialize to a value greater than the largest value in impact info priority code
    count = 0 # initialize counter for anno_id
    if args.anno_type == "snpEff":
        try:
            effect_strings_str = var.INFO["EFF"]
            effect_strings = effect_strings_str.split(",")
            #print len(effect_strings)
        except KeyError:
            if var.INFO.has_key("SNPEFF_EFFECT"):
                impact_details = snpEff.gatk_effect_details(var.INFO)
            else:
                sys.stderr.write("WARNING: The input VCF has no snpEFF annotations. " \
                                 "Variant impact will be set to unknown\n")
        
        for effect_string in effect_strings:
            count += 1
            eff_pieces = snpEff.eff_search.findall(effect_string)
            for piece in eff_pieces:
                impact_string  = piece[0] # the predicted inpact, which is outside the ()
                impact_detail  = piece[1] # all the other information, which is inside the ()
                impact_info    = snpEff.effect_map[impact_string]
                #update the impact stored only if a higher severity transcript is encountered
                if impact_info.priority_code < max_severity:
                    impact_type = impact_info
                    impact_details = snpEff.EffectDetails(impact_string, 
                                                          impact_info.priority, 
                                                          impact_detail, 
                                                          count,
                                                          args.version)
                    max_severity = impact_info.priority_code # store the current "winning" severity for the next iteration.
            
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
                # impact_strings will be [nc_transcript_variant, intron_variant]
                for impact_string in impact_strings:
                    count += 1
                    impact_info = vep.effect_map[impact_string]
                    #update the impact stored only if a higher severity transcript is encountered
                    if impact_info.priority_code < max_severity:
                        impact_type = impact_info
                        impact_details = vep.EffectDetails(impact_string, impact_info.priority, effect_string, count)
                        max_severity = impact_info.priority_code # store the current "winning" severity for the next iteration.
                    
                        
            # we expect VEP to produce a valid impact label for each_string[0]
            elif "&" not in each_string[0]:
                count += 1
                impact_string = each_string[0]
                impact_info = vep.effect_map.get(impact_string)
                if impact_info.priority_code < max_severity:
                    impact_type = impact_info
                    impact_details = vep.EffectDetails(impact_string, impact_info.priority, effect_string, count)
                    #print impact_details.transcript, impact_details.consequence
                    max_severity = impact_info.priority_code # initialize the max_severity to the former value of priority code

    else:
        #should not get here, as the valid -t options should be handled
        #in main()
        sys.exit("ERROR: Unsupported variant annotation type.\n")
        
    return impact_details
