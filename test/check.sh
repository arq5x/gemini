check() {
	echo -n "" "    $3 ... "
	if diff $1 $2 > _err; then
    	echo ok
	else
    	echo fail
		echo $3 >> failed.tests.txt
	fi
}
export -f check
