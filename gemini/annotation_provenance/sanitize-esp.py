import sys
from collections import OrderedDict

j = 0
for line in sys.stdin:
    if line[0] == "#":
        if line.startswith(("##INFO=<ID=EA_AC,", "##INFO=<ID=AA_AC,", "##INFO=<ID=TAC,")):

            line = line.replace(",Number=.", ",Number=R")
            line = line.replace("AltAlleles,RefAllele", "RefAllele,AltAlleles")
            assert "Number=R" in line
        elif line.startswith("##INFO=<ID=GTS,"):
            line = line.replace(",Number=.", ",Number=G")
            assert "Number=G" in line
        print line,
        continue
    j += 1
    # in body, need to adjust GTS, TAC, AA_AC, EA_AC since ESP stores ref last
    # and Number=R means ref should be first.
    fields = line.split("\t")
    info = OrderedDict((p[0], p[1]) for p in (kv.split("=") for kv in fields[7].split(";")))

    # ignore this for now. need to fix if we use them in the db
    #order = info['GTS'].split(",")
    #A1A1,A1A2,A1A3,A1R,A2A2,A2A3,A2R,A3A3,A3R,RR


    alts = fields[4].split(",")
    for field in ("TAC", "AA_AC", "EA_AC"):
        vals = info[field].split(",")
        assert len(vals) == len(alts) + 1, (fields, j)

        vals = vals[-1:] + vals[:-1]
        assert len(vals) == len(alts) + 1, (fields, j)
        info[field] = ",".join(vals)

    k = len(alts) + 1
    for field in ("GTS", "EA_GTC", "AA_GTC"):

        vals = info[field].split(",")
        if len(vals) != k * (k + 1) / 2:
            # X and Y have incorrect numbers here...
            # need to think about this if we end up using them in db
            # but  the GT fields are currently not used.
            assert fields[0] in "XY"

        #info[field] = ",".join(vals[-1:] + vals[:-1])

    fields[7] = ";".join("%s=%s" % kv for kv in info.items())
    print "\t".join(fields),



"""
##fileformat=VCFv4.0
##INFO=<ID=DBSNP,Number=.,Type=String,Description="dbSNP version which established the rs_id">
##INFO=<ID=EA_AC,Number=.,Type=String,Description="European American Allele Count in the order of AltAlleles,RefAllele. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.">
##INFO=<ID=AA_AC,Number=.,Type=String,Description="African American Allele Count in the order of AltAlleles,RefAllele. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.">
##INFO=<ID=TAC,Number=.,Type=String,Description="Total Allele Count in the order of AltAlleles,RefAllele For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.">
##INFO=<ID=MAF,Number=.,Type=String,Description="Minor Allele Frequency in percent in the order of EA,AA,All">
##INFO=<ID=GTS,Number=.,Type=String,Description="Observed Genotypes. For INDELs, A1, A2, or An refers to the N-th alternate allele while R refers to the reference allele.">
##INFO=<ID=EA_GTC,Number=.,Type=String,Description="European American Genotype Counts in the order of listed GTS">
##INFO=<ID=AA_GTC,Number=.,Type=String,Description="African American Genotype Counts in the order of listed GTS">
##INFO=<ID=GTC,Number=.,Type=String,Description="Total Genotype Counts in the order of listed GTS">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Average Sample Read Depth">
##INFO=<ID=FG,Number=.,Type=String,Description="functionGVS">
##INFO=<ID=GM,Number=.,Type=String,Description="accession">
##INFO=<ID=AA,Number=1,Type=String,Description="chimpAllele">
##INFO=<ID=AAC,Number=.,Type=String,Description="aminoAcidChange">
##INFO=<ID=PP,Number=.,Type=String,Description="proteinPosition">
##INFO=<ID=CDP,Number=.,Type=String,Description="cDNAPosition">
##INFO=<ID=PH,Number=.,Type=String,Description="polyPhen">
##INFO=<ID=CP,Number=1,Type=Float,Description="scorePhastCons">
##INFO=<ID=CG,Number=1,Type=Float,Description="consScoreGERP">
##INFO=<ID=GL,Number=.,Type=String,Description="geneList">
##INFO=<ID=GS,Number=.,Type=String,Description="granthamScore">
##INFO=<ID=CA,Number=.,Type=String,Description="clinicalAssociation">
##INFO=<ID=EXOME_CHIP,Number=.,Type=String,Description="Whether a SNP is on the Illumina HumanExome Chip">
##FILTER=<ID=INDEL5,Description="Nearby 1000 Genomes Pilot Indels within 5bp">
##FILTER=<ID=SVM,Description="Failed SVM-based filter at threshold 0.3. (detailed at http://evs.gs.washington.edu/EVS/HelpSNPSummary.jsp#FilterStatus)">
##INFO=<ID=GWAS_PUBMED,Number=.,Type=String,Description="PubMed records for GWAS hits">
##QueryTarget=1:1-249250621
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	69428	rs140739101	T	G	.	PASS	DBSNP=dbSNP_134;EA_AC=313,6535;AA_AC=14,3808;TAC=327,10343;MAF=4.5707,0.3663,3.0647;GTS=GG,GT,TT;EA_GTC=92,129,3203;AA_GTC=1,12,1898;GTC=93,141,5101;DP=110;GL=OR4F5;CP=1.0;CG=0.9;AA=T;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;GM=NM_001005484.1;FG=missense;AAC=CYS/PHE;PP=113/306;CDP=338;GS=205;PH=probably-damaging
1	69476	rs148502021	T	C	.	PASS	DBSNP=dbSNP_134;EA_AC=2,7020;AA_AC=0,3908;TAC=2,10928;MAF=0.0285,0.0,0.0183;GTS=CC,CT,TT;EA_GTC=1,0,3510;AA_GTC=0,0,1954;GTC=1,0,5464;DP=123;GL=OR4F5;CP=0.6;CG=2.3;AA=T;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;GM=NM_001005484.1;FG=missense;AAC=THR/MET;PP=129/306;CDP=386;GS=81;PH=probably-damaging
"""
