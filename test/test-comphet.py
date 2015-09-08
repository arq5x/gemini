import sys
import os
import atexit
import subprocess
from gemini import tests
from gemini.tests import test_inheritance

TestFamily = test_inheritance.TestFamily

HOM_REF, HET, UNKNOWN, HOM_ALT = range(4)

load_cmd = "gemini load -v {name}.vcf -p {name}.ped --skip-gene-tables --test-mode {name}.db"
update_cmd = """echo "UPDATE  variants set is_exonic = 1;" | sqlite3 {name}.db"""


def run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.PIPE)
    except:
        print >>sys.stderr, cmd
        raise


#test_inheritance.main()

fam1 = TestFamily("""\
#family_id  sample_id   paternal_id maternal_id sex phenotype
1   dad   0   0   1  1
1   mom   0   0   2  1
1   kid   dad   mom   1  2
1   kid2   dad   mom   1  1""")


fam1.family.to_ped(open("fam1.ped", "w"))

atexit.register(os.unlink, "fam1.ped")
atexit.register(os.unlink, "fam1.vcf")
atexit.register(os.unlink, "fam1.db")

for ch, gt_types1, gt_types2 in [
        (True, [HOM_REF, HET, HET, HOM_REF], [HET, HOM_REF, HET, HOM_REF]),
        (True, [HET, HET, HET, HOM_REF], [HET, HOM_REF, HET, HOM_REF]),
        (True, [HET, HET, HET, HOM_REF], [HET, HET, HET, HOM_REF]),
        (False, [HET, HET, HET, HOM_REF], [HET, HOM_REF, HET, HOM_ALT]),
        (False, [HOM_REF, HET, HET, HET], [HET, HOM_REF, HET, HET]), # unaffected kit is also CH
        (False, [HOM_ALT, HET, HET, HOM_REF], [HET, HOM_REF, HET, HOM_REF]), # dad is hom_alt
        (False, [HOM_ALT, HOM_REF, HET, HOM_REF], [HOM_REF, HOM_ALT, HET, HOM_REF]), # dad is hom_alt
         ]:

    vfh = open("fam1.vcf", "w")
    fam1.gt_types = gt_types1
    fam1.to_vcf(vfh)
    fam1.gt_types = gt_types2
    fam1.to_vcf(vfh, header=False)
    vfh.close()


    run(load_cmd.format(name="fam1"))
    run(update_cmd.format(name="fam1"))

    ret = run("gemini comp_hets --max-priority 2 --columns 'chrom,start,end,ref,alt' fam1.db")
    if ch:
        assert len(ret.strip().split("\n")) == 2 + 1, ret
    else:
        assert not ret.strip()
    print ret
    print "OK"


fam1 = TestFamily("""\
#family_id  sample_id   paternal_id maternal_id sex phenotype
1   dad   0   0   1  -9
1   mom   0   0   2  -9
1   kid   dad   mom   1  -9""")


fam1.family.to_ped(open("fam2.ped", "w"))

atexit.register(os.unlink, "fam2.ped")
atexit.register(os.unlink, "fam2.vcf")
atexit.register(os.unlink, "fam2.db")

for ch, gt_types1, gt_types2 in [
        (True, [HOM_REF, HET, HET], [HET, HOM_REF, HET]), # the only case that works.
        (True, [HET, HOM_REF, HET], [HOM_REF, HET, HET]), # the only case that works.
        (False, [HET, HET, HET], [HET, HET, HET]),
        (False, [HET, HET, HET], [HET, HOM_REF, HET]),
        (False, [HOM_ALT, HET, HET], [HET, HOM_REF, HET]), # dad is hom_alt
        (False, [HOM_ALT, HOM_REF, HET], [HOM_REF, HOM_ALT, HET]), # dad is hom_alt
         ]:

    vfh = open("fam2.vcf", "w")
    fam1.gt_types = gt_types1
    fam1.to_vcf(vfh)
    fam1.gt_types = gt_types2
    fam1.to_vcf(vfh, header=False)
    vfh.close()


    run(load_cmd.format(name="fam2"))
    run(update_cmd.format(name="fam2"))

    ret = run("gemini comp_hets --pattern-only --columns 'chrom,start,end,ref,alt' fam2.db")
    if ch:
        assert len(ret.strip().split("\n")) == 2 + 1, (ret, ch, gt_types1, gt_types2)
    else:
        assert not ret.strip(), (ret, ch, gt_types1, gt_types2)
    print ret
    print "OK"


famu = TestFamily("""\
#family_id  sample_id   paternal_id maternal_id sex phenotype
1   dad   0   0   1  -9
1   mom   0   0   2  -9
1   kid   dad   mom   1  -9""")
famu.family.to_ped(open("famu.ped", "w"))

vfh = open('famu.vcf', 'w')
famu.gt_types = [UNKNOWN, HET, HET]
famu.to_vcf(vfh)
famu.gt_types = [HET, HOM_REF, HET] # the only case that works.
famu.to_vcf(vfh, header=False)
vfh.close()

run(load_cmd.format(name="famu"))
run(update_cmd.format(name="famu"))

ret = run("gemini comp_hets --pattern-only --columns 'chrom,start,end,ref,alt' famu.db")
assert not len(ret.strip())
ret = run("gemini comp_hets  --columns 'chrom,start,end,ref,alt' famu.db")
assert not len(ret.strip())
