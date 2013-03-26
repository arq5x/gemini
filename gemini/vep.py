#!/usr/bin/env python

#############

# CSQ: Consequence|Codons|Amino_acids|Gene|hgnc|Feature|EXON|polyphen|sift
# non_synonymous_codon|gaT/gaG|D/E|ENSG00000116254|CHD5|ENST00000378006|18/25|benign(0.011)|tolerated(0.3)
# nc_transcript_variant|||ENSG00000116254|CHD5|ENST00000491020|5/6|||
#############

import re
from collections import namedtuple
from collections import defaultdict


class EffectDetails(object):
    def __init__(self, impact_string, severity, detail_string, counter):
        fields = detail_string.split("|")
        self.effect_severity = severity
        self.effect_name = impact_string
        self.anno_id = counter
        self.codon_change = fields[1] if fields[1] != '' else None
        self.aa_change = fields[2] if fields[2] != '' else None
        self.ensembl_gene = fields[3] if fields[3] != '' else None
        self.hgnc = fields[4] if fields[4] != '' else None
        self.gene = self.hgnc if fields[4] != '' else self.ensembl_gene
        self.transcript = fields[5] if fields[5] != '' else None
        self.exon = fields[6] if fields[6] != '' else None
        self.polyphen = fields[7] if fields[7] != '' else None
        self.sift = fields[8] if fields[8] != '' else None
        self.aa_length = None
        self.biotype = None
        self.warnings = None
        self.consequence = effect_dict[
            self.effect_name] if self.effect_severity != None else self.effect_name
        if len(fields) > 9:
            self.warnings = fields[9]

        # rules for being exonic.
        # 1. the impact must be in the list of exonic impacts
        # 2. the exon information != None
        self.is_exonic = 0
        if self.effect_name in exonic_impacts and \
                self.exon is not None:
            self.is_exonic = 1

        self.is_lof = 0 if self.effect_severity != "HIGH" else 1
        # Exons that are coding (excludes UTR's)
        if self.is_exonic == 0:
            self.is_coding = 0
        elif self.is_exonic == 1:
            if self.effect_name.startswith("5_prime_UTR") or self.effect_name.startswith("3_prime_UTR"):
                self.is_coding = 0
            else:
                self.is_coding = 1
        # parse Polyphen predictions
        if self.polyphen is not None:
            self.polyphen_b = self.polyphen.split("(")
            self.polyphen_pred = self.polyphen_b[0]
            self.polyphen2 = self.polyphen_b[1].split(")")
            self.polyphen_score = self.polyphen2[0]
        else:
            self.polyphen_pred = None
            self.polyphen_score = None
        # parse SIFT predictions
        if self.sift is not None:
            self.sift_b = self.sift.split("(")
            self.sift_pred = self.sift_b[0]
            self.sift2 = self.sift_b[1].split(")")
            self.sift_score = self.sift2[0]
        else:
            self.sift_pred = None
            self.sift_score = None

    def __str__(self):

        return "\t".join([self.consequence, self.effect_severity, self.codon_change,
                          self.aa_change, self.aa_length, self.biotype,
                          self.ensembl_gene, self.gene, self.transcript,
                          self.exon, self.is_exonic, self.anno_id, self.polyphen_pred,
                          self.polyphen_score, self.sift_pred, self.sift_score,
                          self.is_coding, self.is_lof])

    def __repr__(self):
        return self.__str__()


exonic_impacts = ["stop_gained",
                  "stop_lost",
                  "frameshift_variant",
                  "initiator_codon_variant",
                  "inframe_deletion",
                  "inframe_insertion",
                  "missense_variant",
                  "incomplete_terminal_codon_variant",
                  "stop_retained_variant",
                  "synonymous_variant",
                  "coding_sequence_variant",
                  "5_prime_UTR_variant",
                  "3_prime_UTR_variant",
                  "transcript_ablation",
                  "transcript_amplification",
                  "feature_elongation",
                  "feature_truncation"]


effect_names = ["splice_acceptor_variant", "splice_donor_variant",
                "stop_gained", "stop_lost",
                "non_coding_exon_variant", "frameshift_variant",
                "initiator_codon_variant", "inframe_deletion",
                "inframe_insertion", "missense_variant",
                "splice_region_variant", "incomplete_terminal_codon_variant",
                "stop_retained_variant", "synonymous_variant",
                "coding_sequence_variant", "mature_miRNA_variant",
                "5_prime_UTR_variant", "3_prime_UTR_variant",
                "intron_variant", "NMD_transcript_variant",
                "nc_transcript_variant", "upstream_gene_variant",
                "downstream_gene_variant", "regulatory_region_variant",
                "TF_binding_site_variant", "intergenic_variant",
                "regulatory_region_ablation", "regulatory_region_amplification",
                "transcript_ablation", "transcript_amplification",
                "TFBS_ablation", "TFBS_amplification",
                "feature_elongation", "feature_truncation"]


effect_dict = defaultdict()
effect_dict = {
    'splice_acceptor_variant': 'splice_acceptor', 'splice_donor_variant': 'splice_donor',
    'stop_gained': 'stop_gain', 'stop_lost': 'stop_loss',
    'non_coding_exon_variant': 'nc_exon', 'frameshift_variant': 'frame_shift',
    'initiator_codon_variant': 'transcript_codon_change', 'inframe_deletion': 'inframe_codon_loss',
    'inframe_insertion': 'inframe_codon_gain', 'missense_variant': 'non_syn_coding',
    'splice_region_variant': 'other_splice_variant', 'incomplete_terminal_codon_variant': 'incomplete_terminal_codon',
    'stop_retained_variant': 'synonymous_stop', 'synonymous_variant': 'synonymous_coding',
    'coding_sequence_variant': 'CDS', 'mature_miRNA_variant': 'mature_miRNA',
    '5_prime_UTR_variant': 'UTR_5_prime', '3_prime_UTR_variant': 'UTR_3_prime',
    'intron_variant': 'intron', 'NMD_transcript_variant': 'NMD_transcript',
    'nc_transcript_variant': 'nc_transcript', 'upstream_gene_variant': 'upstream',
    'downstream_gene_variant': 'downstream', 'regulatory_region_variant': 'regulatory_region',
    'TF_binding_site_variant': 'TF_binding_site', 'intergenic_variant': 'intergenic',
    'regulatory_region_ablation': 'regulatory_region_ablation', 'regulatory_region_amplification': 'regulatory_region_amplification',
    'transcript_ablation': 'transcript_ablation', 'transcript_amplification': 'transcript_amplification',
    'TFBS_ablation': 'TFBS_ablation', 'TFBS_amplification': 'TFBS_amplification',
    'feature_elongation': 'feature_elongation', 'feature_truncation': 'feature_truncation'}

effect_desc = ["The variant hits the splice acceptor site (2 basepair region at 3' end of an intron)", "The variant hits the splice donor site (2 basepair region at 5'end of an intron)",
               "Variant causes a STOP codon", "Variant causes stop codon to be mutated into a non-stop codon",
               "Variant causes a change in the non coding exon sequence", "Insertion or deletion causes a frame shift in coding sequence",
               "Variant causes atleast one base change in the first codon of a transcript", "An inframe non-syn variant that deletes bases from the coding sequence",
               "An inframe non-syn variant that inserts bases in the coding sequence", "The variant causes a different amino acid in the coding sequence",
               "Variant causes a change within the region of a splice site (1-3bps into an exon or 3-8bps into an intron)", "The variant hits the incomplete codon of a transcript whose end co-ordinate is not known",
               "The variant causes stop codon to be mutated into another stop codon", "The variant causes no amino acid change in coding sequence",
               "Variant hits coding sequence with indeterminate effect", "The variant hits a microRNA",
               "Variant hits the 5 prime untranslated region", "Variant hits the 3 prime untranslated region",
               "Variant hits an intron", "A variant hits a transcript that is predicted to undergo nonsense mediated decay",
                   "Variant hits a gene that does not code for a protein", "The variant hits upstream of a gene (5' of a gene)",
                   "The variant hits downstream of a gene (3' of a gene)", "Variant hits the regulatory region annotated by Ensembl(e.g promoter)",
                   "Variant falls in a transcription factor binding motif within an Ensembl regulatory region", "The variant is located in the intergenic region, between genes",
                   "SV causes ablation of a regulatory region", "SV results in an amplification of a regulatory region",
                   "SV causes an ablation/deletion of a transcript feature", "SV causes an amplification of a transcript feature",
                   "SV results in a deletion of the TFBS", "SV results in an amplification of a region containing TFBS",
                   "SV causes an extension of a genomic feature wrt reference", "SV causes a reduction of a genomic feature compared to reference"]

effect_priorities = ["HIGH", "HIGH",
                     "HIGH", "HIGH",
                     "LOW", "HIGH",
                     "HIGH", "MED",
                     "MED", "MED",
                     "MED", "LOW",
                     "LOW", "LOW",
                     "LOW", "MED",
                     "LOW", "LOW",
                     "LOW", "LOW",
                     "LOW", "LOW",
                     "LOW", "MED",
                     "MED", "LOW",
                     "MED", "MED",
                     "LOW", "LOW",
                     "MED", "MED",
                     "LOW", "LOW"]

effect_priority_codes = [1, 1,
                        1, 1,
                        3, 1,
                        1, 2,
                        2, 2,
                        2, 3,
                        3, 3,
                        3, 2,
                        3, 3,
                        3, 3,
                        3, 3,
                        3, 2,
                        2, 3,
                        2, 2,
                        3, 3,
                        2, 2,
                        3, 3]

effect_ids = range(1, 35)

effect_map = {}
EffectInfo = namedtuple(
    'EffectInfo', ['id', 'priority', 'priority_code', 'desc'])
for i, effect_name in enumerate(effect_names):
    info = EffectInfo(effect_ids[i], effect_priorities[i],
                      effect_priority_codes[i], effect_desc[i])
    effect_map[effect_name] = info
