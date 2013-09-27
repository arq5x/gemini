#!/usr/bin/env python
#############################################
# consequence as per snpeff v3:
# Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon
# NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aCg/aTg|T143M|459|XKR3|protein_coding|CODING|ENST00000331428|exon_22_17280661_17280914

#############################################

import re
from collections import namedtuple
from collections import defaultdict


class EffectDetails(object):
    def __init__(self, name, severity, detail_string, counter, snp_eff_version):
        fields = detail_string.split("|")
        self.effect_name = name
        self.anno_id = counter
        self.effect_severity = severity
        self.impact = fields[1] if fields[1] != '' else None
        self.codon = fields[2] if fields[2] != '' else None
        self.aa_change = fields[3] if fields[3] != '' else None

        # snpEff >= v3.0 includes aa_length
        if snp_eff_version is not None and snp_eff_version >= 3:
            self.aa_length = fields[4] if fields[4] != '' else None
            self.gene = fields[5] if fields[5] != '' else None
            self.biotype = fields[6] if fields[6] != '' else None
            self.coding = fields[7] if fields[7] != '' else None
            self.transcript = fields[8] if fields[8] != '' else None
            self.exon = fields[9] if fields[9] != '' else None
            self.warnings = None
            if len(fields) > 9:
                self.warnings = fields[9]
        else:
            self.aa_length = None
            self.gene = fields[4] if fields[4] != '' else None
            self.biotype = fields[5] if fields[5] != '' else None
            self.coding = fields[6] if fields[6] != '' else None
            self.transcript = fields[7] if fields[7] != '' else None
            self.exon = fields[8] if fields[8] != '' else None
            self.warnings = None
            if len(fields) > 8:
                self.warnings = fields[8]

        # Handling only codon change and not distance(snpEff v3.3)
        if self.effect_name not in ("DOWNSTREAM", "UPSTREAM"):
            self.codon_change = self.codon
        else:
            self.codon_change =  None
        
        # rules for being exonic.
        # 1. must be protein_coding
        # 2. the impact must be in a list of impacts
        #    that are known to be exonic.
        if self.biotype == "protein_coding" and \
                self.effect_name in exonic_impacts:
            self.is_exonic = 1
        else:
            self.is_exonic = 0

        # rules for being exonic.
        # 1. must be protein_coding
        # 2. must be exonic, yet must not be a UTR
        self.is_coding = 0
        if self.is_exonic and not (self.effect_name == "START_GAINED" or
                                   self.effect_name.startswith("UTR_")):
            self.is_coding = 1

        # rules for being loss-of-function (lof).
        # 1. must be protein_coding
        # 2. must be a coding variant with HIGH impact
        # 3. must affect a protein_coding transcript
        self.is_lof = 0
        if self.effect_severity == "HIGH" and self.biotype == "protein_coding":
            self.is_lof = 1

        self.polyphen_pred = None
        self.polyphen_score = None
        self.sift_pred = None
        self.sift_score = None
        self.consequence = effect_dict[self.effect_name]

    def __str__(self):
        return "\t".join([self.consequence, self.effect_severity,
                          str(self.impact), str(self.codon_change),
                          str(self.aa_change), str(self.aa_length), str(self.gene),
                          str(self.biotype), str(self.is_exonic),
                          str(self.is_coding), str(self.transcript),
                          str(self.exon), str(self.anno_id)])

    def __repr__(self):
        return self.__str__()


exonic_impacts = ["CODON_CHANGE",
                  "CODON_CHANGE_PLUS_CODON_DELETION",
                  "CODON_CHANGE_PLUS_CODON_INSERTION",
                  "CODON_DELETION",
                  "CODON_INSERTION",
                  "EXON",
                  "EXON_DELETED",
                  "FRAME_SHIFT",
                  "GENE",
                  "NON_SYNONYMOUS_CODING",
                  "RARE_AMINO_ACID",
                  "START_GAINED",
                  "START_LOST",
                  "STOP_GAINED",
                  "STOP_LOST",
                  "SYNONYMOUS_CODING",
                  "SYNONYMOUS_START",
                  "SYNONYMOUS_STOP",
                  "TRANSCRIPT",
                  "UTR_3_DELETED",
                  "UTR_3_PRIME",
                  "UTR_5_DELETED",
                  "UTR_5_PRIME",
                  "NON_SYNONYMOUS_START"]


effect_names = ["CDS",
                "CODON_CHANGE",
                "CODON_CHANGE_PLUS_CODON_DELETION",
                "CODON_CHANGE_PLUS_CODON_INSERTION",
                "CODON_DELETION",
                "CODON_INSERTION",
                "DOWNSTREAM",
                "EXON",
                "EXON_DELETED",
                "FRAME_SHIFT",
                "GENE",
                "INTERGENIC",
                "INTERGENIC_CONSERVED",
                "INTRAGENIC",
                "INTRON",
                "INTRON_CONSERVED",
                "NON_SYNONYMOUS_CODING",
                "RARE_AMINO_ACID",
                "SPLICE_SITE_ACCEPTOR",
                "SPLICE_SITE_DONOR",
                "START_GAINED",
                "START_LOST",
                "STOP_GAINED",
                "STOP_LOST",
                "SYNONYMOUS_CODING",
                "SYNONYMOUS_START",
                "SYNONYMOUS_STOP",
                "TRANSCRIPT",
                "UPSTREAM",
                "UTR_3_DELETED",
                "UTR_3_PRIME",
                "UTR_5_DELETED",
                "UTR_5_PRIME",
                "NON_SYNONYMOUS_START",
                "NONE"]

effect_dict = defaultdict()
effect_dict = {'CDS': 'CDS',
               'CODON_CHANGE': 'inframe_codon_change',
               'CODON_CHANGE_PLUS_CODON_DELETION': 'codon_change_del',
               'CODON_CHANGE_PLUS_CODON_INSERTION': 'codon_change_ins',
               'CODON_DELETION': 'inframe_codon_loss',
               'CODON_INSERTION': 'inframe_codon_gain',
               'DOWNSTREAM': 'downstream',
               'EXON': 'exon',
               'EXON_DELETED': 'exon_deleted',
               'FRAME_SHIFT': 'frame_shift',
               'GENE': 'gene',
               'INTERGENIC': 'intergenic',
               'INTERGENIC_CONSERVED': 'intergenic',
               'INTRAGENIC': 'intragenic',
               'INTRON': 'intron',
               'INTRON_CONSERVED': 'intron_conserved',
               'NON_SYNONYMOUS_CODING': 'non_syn_coding',
               'RARE_AMINO_ACID': 'rare_amino_acid',
               'SPLICE_SITE_ACCEPTOR': 'splice_acceptor',
               'SPLICE_SITE_DONOR': 'splice_donor',
               'START_GAINED': 'start_gain',
               'START_LOST': 'start_loss',
               'STOP_GAINED': 'stop_gain',
               'STOP_LOST': 'stop_loss',
               'SYNONYMOUS_CODING': 'synonymous_coding',
               'SYNONYMOUS_START': 'synonymous_start',
               'SYNONYMOUS_STOP': 'synonymous_stop',
               'TRANSCRIPT': 'transcript',
               'UPSTREAM': 'upstream',
               'UTR_3_DELETED': 'UTR_3_del',
               'UTR_3_PRIME': 'UTR_3_prime',
               'UTR_5_DELETED': 'UTR_5_del',
               'UTR_5_PRIME': 'UTR_5_prime',
               'NON_SYNONYMOUS_START': 'non_synonymous_start',
               'NONE': 'None'}

effect_desc = ["The variant hits a CDS.",
               "One or many codons are changed",
               "One codon is changed and one or more codons are deleted",
               "One codon is changed and one or many codons are inserted",
               "One or many codons are deleted",
               "One or many codons are inserted",
               "Downstream of a gene (default length: 5K bases)",
               "The variant hits an exon.",
               "A deletion removes the whole exon.",
               "Insertion or deletion causes a frame shift",
               "The variant hits a gene.",
               "The variant is in an intergenic region.",
               "The variant is in a highly conserved intergenic region.",
               "The variant hits a gene, but no transcripts within \
               the gene.",
               "Variant hits an intron.",
               "The variant is in a highly conserved intronic region.",
               "Variant causes a codon that produces a different \
               amino acid.",
               "The variant hits a rare amino acid thus is likely to \
                produce protein loss of function",
               "The variant hits a splice acceptor site (defined as two \
               bases before exon start, except for the first exon).",
               "The variant hits a Splice donor site (defined as two \
               bases after coding exon end, except for the last exon).",
               "A variant in 5'UTR region produces a three base sequence \
                that can be a START codon.",
               "Variant causes start codon to be mutated into a non-start \
                codon.",
               "Variant causes a STOP codon.",
               "Variant causes stop codon to be mutated into a non-stop \
                codon.",
               "Variant causes a codon that produces the same amino acid",
               "Variant causes start codon to be mutated into another \
               start codon.",
               "Variant causes stop codon to be mutated into another stop \
                codon.",
               "The variant hits a transcript.",
               "Upstream of a gene (default length: 5K bases).",
               "The variant deletes and exon which is in the 3'UTR of the \
                transcript.",
               "Variant hits 3'UTR region.",
               "The variant deletes and exon which is in the 5'UTR of the \
                transcript.",
               "Variant hits 5'UTR region.",
               "The variant causes a start codon to be changed into a \
                different codon",
                "Unknown"]

effect_priorities = ["LOW",
                     "MED",
                     "MED",
                     "MED",
                     "MED",
                     "MED",
                     "LOW",
                     "LOW",
                     "HIGH",
                     "HIGH",
                     "LOW",
                     "LOW",
                     "LOW",
                     "LOW",
                     "LOW",
                     "LOW",
                     "MED",
                     "HIGH",
                     "HIGH",
                     "HIGH",
                     "LOW",
                     "HIGH",
                     "HIGH",
                     "HIGH",
                     "LOW",
                     "LOW",
                     "LOW",
                     "LOW",
                     "LOW",
                     "MED",
                     "LOW",
                     "MED",
                     "LOW",
                     "HIGH",
                     "LOW"]

effect_priority_codes = [3,
                         2,
                         2,
                         2,
                         2,
                         2,
                         3,
                         3,
                         1,
                         1,
                         3,
                         3,
                         3,
                         3,
                         3,
                         3,
                         2,
                         1,
                         1,
                         1,
                         3,
                         1,
                         1,
                         1,
                         3,
                         3,
                         3,
                         3,
                         3,
                         2,
                         3,
                         2,
                         3,
                         1,
                         3]

effect_ids = range(1, 36)
effect_map = {}
EffectInfo = namedtuple(
    'EffectInfo', ['id', 'priority', 'priority_code', 'desc'])

for i, effect_name in enumerate(effect_names):
    info = EffectInfo(effect_ids[i], effect_priorities[i],
                      effect_priority_codes[i], effect_desc[i])
    effect_map[effect_name] = info

eff_pattern = '(\S+)[(](\S+)[)]'
eff_search = re.compile(eff_pattern)


def gatk_effect_details(info):
    """Convert GATK prepared snpEff effect details into standard EffectDetails.
    """
    name = info.get("SNPEFF_EFFECT", None)
    if name is not None:
        effect = effect_map[name]
        detail_string = "|{impact}|{codon_change}|{aa_change}|{gene}|{biotype}|{coding}|{transcript}|{exon}".format(
            impact=info.get("SNPEFF_IMPACT", ""),
            codon_change=info.get("SNPEFF_CODON_CHANGE", ""),
            aa_change=info.get("SNPEFF_AMINO_ACID_CHANGE", ""),
            gene=info.get("SNPEFF_GENE_NAME", ""),
            biotype=info.get("SNPEFF_GENE_BIOTYPE", ""),
            coding="",
            transcript=info.get("SNPEFF_TRANSCRIPT", ""),
            exon=info.get("SNPEFF_EXON_ID", ""))
        return EffectDetails(name, effect.priority, detail_string, 0)
