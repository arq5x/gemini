check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
        exit 1
    fi
}

echo "test -t.all.1"
# annotated with only snpEff, make sure -t all doesn't fail
gemini load -t all -v test.snpeff.vcf --test-mode --skip-gene-tables test.tall.snpeff.db 2> err
echo "FAM138A
FAM138A
FAM138A
FAM138A
FAM138A
OR4F5
OR4F5
OR4F5
OR4F5
OR4F5" > exp
gemini query -q "select gene from variants" test.tall.snpeff.db  > obs
check obs exp
rm obs exp

# annotated with only VEP, make sure -t all doesn't fail
echo "test -t.all.2"
gemini load -t all -v  test-vep-extra.vcf --test-mode --skip-gene-tables test.tall.vep.db 2> err
echo "NFE2L2
None" > exp
gemini query -q "select gene from variants" test.tall.vep.db > obs
check obs exp
rm obs exp
