#!/usr/bin/env python

#############

#CSQ: Consequence|Codons|Amino_acids|Gene|hgnc|Feature|EXON|polyphen|sift|condel
#non_synonymous_codon|gaT/gaG|D/E|ENSG00000116254|CHD5|ENST00000378006|18/25|benign(0.011)|tolerated(0.3)|neutral(0.029)
#nc_transcript_variant|||ENSG00000116254|CHD5|ENST00000491020|5/6|||

#############

import re
from collections import namedtuple
from collections import defaultdict

class EffectDetails(object):
    def __init__(self, impact_string, severity, detail_string):
        fields = detail_string.split("|")
        self.effect_severity = severity
        self.effect_name = impact_string
        self.codon_change = fields[1] if fields[1] != '' else None
        self.aa_change = fields[2] if fields[2] != '' else None
        self.ensembl_gene = fields[3] if fields[3] != '' else None
        self.gene = fields[4] if fields[4] != '' else None
        self.transcript = fields[5] if fields[5] != '' else None
        self.exon = fields[6] if fields[6] != '' else None
        self.polyphen =  fields[7] if fields[7] != '' else None
        self.sift = fields[8] if fields[8] != '' else None
        self.condel = fields[9] if fields[9] != '' else None
        self.warnings = None
        self.consequence = effect_dict[self.effect_name]
        if len(fields) > 9:
            self.warnings = fields[9]
        
        self.exonic = 0 if self.exon is None else 1
        self.is_lof = 0 if self.effect_severity != "HIGH" else 1
        # Exons that are coding (excludes UTR's)
        if self.exonic == 0:
            self.coding = 0 
        elif self.exonic == 1:
            if self.effect_name.startswith("5_prime_UTR") or self.effect_name.startswith("3_prime_UTR"):
                self.coding = 0 
            else:
                self.coding = 1

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
        # parse CONDEL predictions
        if self.condel is not None:
                self.condel_b = self.condel.split("(")
                self.condel_pred = self.condel_b[0]
                self.condel2 = self.condel_b[1].split(")")
                self.condel_score = self.condel2[0]
        else:
            self.condel_pred  = None
            self.condel_score = None

    def __str__(self):

        return "\t".join([self.consequence, self.effect_severity,
                          self.codon_change, self.aa_change,
                          self.ensembl_gene, self.gene, self.transcript,
                          self.exon, self.exonic, self.polyphen_pred,
                          self.polyphen_score, self.sift_pred, self.sift_score,
                          self.condel_pred, self.condel_score, self.coding, self.is_lof])
    def __repr__(self):
        return self.__str__()



effect_names    = ["splice_acceptor_variant", "splice_donor_variant", 
                   "stop_gained", "stop_lost", 
                   "complex_change_in_transcript", "frameshift_variant",
                   "initiator_codon_change", "inframe_codon_loss",
                   "inframe_codon_gain","non_synonymous_codon",
                   "splice_region_variant", "incomplete_terminal_codon_variant",
                   "stop_retained_variant", "synonymous_codon",
                   "coding_sequence_variant", "mature_miRNA_variant",
                   "5_prime_UTR_variant", "3_prime_UTR_variant",
                   "intron_variant", "NMD_transcript_variant",
                   "nc_transcript_variant", "2KB_upstream_variant",
                   "5KB_upstream_variant", "500B_downstream_variant",
                   "5KB_downstream_variant", "regulatory_region_variant",
                   "TF_binding_site_variant", "intergenic_variant",
                   "entirely_within_Gene", "deletion", "entirely_within_Transcript", "partial_overlap"]
                   

effect_dict = defaultdict()
effect_dict = {'splice_acceptor_variant': 'splice_acceptor', 'splice_donor_variant': 'splice_donor', 'stop_gained': 'stop_gain', 
              'stop_lost': 'stop_loss', 'complex_change_in_transcript': 'complex_change-transcript', 'frameshift_variant': 'frame_shift', 
              'initiator_codon_change': 'transcript_codon_change', 'inframe_codon_loss': 'inframe_codon_loss', 'inframe_codon_gain': 'inframe_codon_gain',
              'non_synonymous_codon': 'non_syn_coding', 'splice_region_variant': 'other_splice_variant', 'incomplete_terminal_codon_variant': 'incomplete_terminal_codon',
              'stop_retained_variant': 'synonymous_stop', 'synonymous_codon': 'synonymous_coding', 'coding_sequence_variant': 'CDS',
              'mature_miRNA_variant': 'mature_miRNA', '5_prime_UTR_variant': 'UTR_5_prime', '3_prime_UTR_variant': 'UTR_3_prime', 
              'intron_variant': 'intron', 'NMD_transcript_variant': 'NMD_transcript', 'nc_transcript_variant': 'nc_transcript',
              '2KB_upstream_variant': 'upstream', '5KB_upstream_variant': 'upstream', '500B_downstream_variant': 'downstream', 
              '5KB_downstream_variant': 'downstream', 'regulatory_region_variant': 'regulatory_region', 'TF_binding_site_variant': 'TF_binding_site',
              'intergenic_variant': 'intergenic', 'entirely_within_Gene': 'entirely_within_Gene', 'deletion': 'deletion',
              'entirely_within_Transcript': 'entirely_within_Transcript', 'partial_overlap': 'partial_overlap' };

effect_desc     = ["The variant hits the splice acceptor site (2 basepair region at 3' end of an intron)", "The variant hits the splice donor site (2 basepair region at 5'end of an intron)",
                   "Variant causes a STOP codon", "Variant causes stop codon to be mutated into a non-stop codon", 
                   "Variant causes an Indel spanning a exon/intron or coding/UTR border", "Insertion or deletion causes a frame shift in coding sequence", 
                   "Variant causes atleast one base change in the first codon of a transcript", "The variant causes an inframe deletion of a codon",
                   "The variant causes an inframe insertion of a new codon", "The variant causes a different amino acid in the coding sequence",
                   "Variant causes a change within the region of a splice site (1-3bps into an exon or 3-8bps into an intron)", "The variant hits the incomplete codon of a transcript whose end co-ordinate is not known",
                   "The variant causes stop codon to be mutated into another stop codon", "The variant causes no amino acid change in coding sequence",
                   "Variant hits coding sequence with indeterminate effect", "The variant hits a microRNA",
                   "Variant hits the 5 prime untranslated region", "Variant hits the 3 prime untranslated region",
                   "Variant hits an intron", "A variant hits a transcript that is predicted to undergo nonsense mediated decay",
                   "Variant hits a gene that does not code for a protein", "The variant hits upstream of a gene (within 2KB 5' of a gene)",
                   "Variant hits upstream of a gene (within 5KB 5' of a gene)", "The variant hits downstream of a gene (within half kb of the end of gene)", 
                   "The variant hits downstream of a gene (within 5KB of the end of gene)", "Variant hits the regulatory region annotated by Ensembl(e.g promoter)", 
                   "Variant falls in a transcription factor binding motif within an Ensembl regulatory region", "The variant is located in the intergenic region, between genes",
                   "SV within gene", "SV- a deletion", "SV within transcript", "SV - partial overlap"]

effect_priorities = ["HIGH", "HIGH",
                   "HIGH", "HIGH",
                   "MED", "HIGH",
                   "HIGH", "MED",
                   "MED", "MED",
                   "MED", "LOW",
                   "LOW", "LOW",
                   "LOW", "MED", 
                   "LOW", "LOW",
                   "LOW", "LOW",
                   "LOW", "LOW",
                   "LOW", "LOW",
                   "LOW", "MED",
                   "MED", "LOW",
                   "LOW", "LOW",
                   "MED", "LOW"]

effect_priority_codes = [1, 1,
                        1, 1,
                        2, 1,
                        1, 2,
                        2, 2,
                        2, 3,
                        3, 3,
                        3, 2, 
                        3, 3,
                        3, 3,
                        3, 3,
                        3, 3,
                        3, 2,
                        2, 3,
                        3, 3,
                        2, 3]

effect_ids = range(1,33)

effect_map = {}
EffectInfo = namedtuple('EffectInfo', ['id', 'priority', 'priority_code', 'desc'])
for i, effect_name in enumerate(effect_names):
    info = EffectInfo(effect_ids[i], effect_priorities[i], effect_priority_codes[i], effect_desc[i])
    effect_map[effect_name] = info