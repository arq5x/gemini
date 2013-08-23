gemini load -v test.snpeff.vcf -t snpEff test.snpeff.vcf.db
gemini load -v test1.snpeff.vcf -t snpEff test1.snpeff.db
gemini load -v test2.snpeff.vcf test2.snpeff.db
gemini load -v test3.snpeff.vcf test3.snpeff.db
gemini load -v test.clinvar.vcf test.clinvar.db
gemini load -v test4.vep.snpeff.vcf -t snpEff test4.snpeff.db
gemini load -v test4.vep.snpeff.vcf -t VEP test4.vep.db
gemini load -v test5.vep.snpeff.vcf -t snpEff test5.snpeff.db
gemini load -v test5.vep.snpeff.vcf -t VEP test5.vep.db
gemini load -v test.query.vcf -t snpEff test.query.db
gemini load -v test.region.vep.vcf -t VEP test.region.db
gemini load -v test.burden.vcf -t VEP -p test.burden.ped test.burden.db
gemini load -v test.auto_dom.vcf -t snpEff -p test.auto_dom.ped test.auto_dom.db
gemini load -v test.auto_rec.vcf -t snpEff -p test.auto_rec.ped test.auto_rec.db
