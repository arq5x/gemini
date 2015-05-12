for f in *.db; do
	echo $f
	if [[ $f != *1000G* ]]; then
	   gemini bcolz_index $f;
	fi
done
