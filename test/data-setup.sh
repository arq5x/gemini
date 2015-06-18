gemini load -v test.auto_rec.vcf tt.db --no-genotypes && rm tt.db
gemini load --skip-gene-tables --test-mode --skip-gerp --skip-cadd -v test.PLs.vcf  test.PLs.db --no-bcolz
gemini load --skip-gene-tables --test-mode --skip-gerp --skip-cadd -v test.mendel.vcf -p test.mendel.ped  test.mendel.db --no-bcolz
gemini load --skip-gene-tables --test-mode --skip-gerp --skip-cadd -v test.multiple-alts.decomp.snpeff.vcf -t snpEff test.mad.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.cosmic.vcf --skip-gerp-bp --skip-cadd test.cosmic.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.roh.vcf --skip-gerp-bp --skip-cadd -t snpEff test.roh.vcf.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.cadd.vcf --skip-gerp-bp --skip-cadd test.no.cadd.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.snpeff.vcf --skip-gerp-bp --skip-cadd -t snpEff test.snpeff.vcf.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test1.snpeff.vcf --skip-gerp-bp --skip-cadd -t snpEff test1.snpeff.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test1.snpeff.vcf --skip-gerp-bp --skip-cadd -t VEP test1.vep.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test2.snpeff.vcf --skip-gerp-bp --skip-cadd test2.snpeff.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test3.snpeff.vcf --skip-gerp-bp --skip-cadd test3.snpeff.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.clinvar.vcf --skip-gerp-bp --skip-cadd test.clinvar.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test4.vep.snpeff.vcf --skip-gerp-bp --skip-cadd -t snpEff test4.snpeff.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test4.vep.snpeff.vcf --skip-gerp-bp --skip-cadd -t VEP test4.vep.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test5.vep.snpeff.vcf --skip-gerp-bp --skip-cadd -t snpEff test5.snpeff.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test5.vep.snpeff.vcf --skip-gerp-bp --skip-cadd -t VEP test5.vep.db --no-bcolz
gemini load --test-mode -v test.query.vcf --skip-gerp-bp --skip-cadd -t snpEff test.query.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.query.vcf --skip-gerp-bp --skip-cadd -t VEP test.query.vep.db --no-bcolz
gemini load --test-mode -v test.query.vcf --skip-gerp-bp --skip-cadd --cores 2 -t snpEff test.query.core.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.region.vep.vcf --skip-gerp-bp --skip-cadd -t VEP test.region.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.burden.vcf --skip-gerp-bp --skip-cadd -t VEP -p test.burden.ped test.burden.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.ped test.comp_het.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.2.vcf --skip-gerp-bp --skip-cadd -t snpEff test.comp_het_default.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.3.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.ped test.comp_het_default.2.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.4.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.2.ped test.comp_het_default.4.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.5.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.2.ped test.comp_het_default.5.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.6.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.2.ped test.comp_het_default.6.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.comp_het.7.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.2.ped test.comp_het_default.7.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.ped test.auto_dom.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.no_parents.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.no_parents.ped test.auto_dom.no_parents.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.no_parents.2.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.no_parents.2.ped test.auto_dom.no_parents.2.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.no_parents.3.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.no_parents.3.ped test.auto_dom.no_parents.3.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.no_parents.4.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.no_parents.4.ped test.auto_dom.no_parents.4.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_dom.no_parents.5.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_dom.no_parents.5.ped test.auto_dom.no_parents.5.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.ped test.auto_rec.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.no_parents.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.no_parents.ped test.auto_rec.no_parents.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.no_parents.2.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.no_parents.2.ped test.auto_rec.no_parents.2.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.no_parents.3.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.no_parents.3.ped test.auto_rec.no_parents.3.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.no_parents.4.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.no_parents.4.ped test.auto_rec.no_parents.4.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.auto_rec.no_parents.5.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.auto_rec.no_parents.5.ped test.auto_rec.no_parents.5.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.de_novo.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.de_novo.ped test.de_novo.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.de_novo.affected.and.unaffected.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.de_novo.affected.and.unaffected.ped test.de_novo.affected.and.unaffected.db --no-bcolz
gemini load --skip-gene-tables --test-mode -p test4.snpeff.ped -v test4.vep.snpeff.vcf --skip-gerp-bp --skip-cadd -t snpEff test4.snpeff.ped.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.vcf_id.snpeff.vcf  --skip-gerp-bp --skip-cadd -t snpEff test.vcf_id.snpeff.vcf.db --no-bcolz
gemini load --skip-gene-tables --test-mode -p test.de_novo.ped -v test.family.vcf  --skip-gerp-bp --skip-cadd -t snpEff test.family.db --no-bcolz
gemini load --skip-gene-tables --test-mode -p test_extended_ped.ped -v test4.vep.snpeff.vcf  --skip-gerp-bp --skip-cadd -t snpEff extended_ped.db --no-bcolz
cp extended_ped.db test.amend.db
gemini load -p test.somatic.ped -v test.somatic.vcf  --skip-gerp-bp --skip-cadd -t VEP test.somatic.db --no-bcolz
gemini load -p test.fusions.ped -v test.fusions.vcf  --skip-gerp-bp --skip-cadd -t VEP test.fusions.db --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.exac.vcf --skip-gerp-bp --skip-cadd test.exac.db --no-bcolz
gemini load -t VEP -v test.esp.zero.vcf test.esp.db --skip-gerp-bp --skip-gene-tables --test-mode --no-bcolz
gemini load --skip-gene-tables --test-mode -v test.dashes.vcf --skip-gerp-bp --skip-cadd -t snpEff test.dashes.db -p test.dashes.ped

# test-phasing
gemini load --skip-gene-tables --test-mode -v test.comp_het.unphase.vcf --skip-gerp-bp --skip-cadd -t snpEff -p test.comp_het.ped test.comp_het.unphase.db --no-bcolz
