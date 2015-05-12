check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

for f in test*.sh; do
	sed 's/gemini query/gemini query --use-bcolz/' $f | bash
done
