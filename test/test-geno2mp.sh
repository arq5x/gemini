check()
{
    if diff $1 $2; then
        echo ok
    else
        echo fail
    fi
}
export -f check

echo "geno2mp.t1..."
gemini query -q "select distinct(geno2mp_hpo_ct) from variants" test.geno2mp.db > obs
echo "5
-1" > exp
check obs exp


