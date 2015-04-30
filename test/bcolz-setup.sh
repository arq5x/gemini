for f in *.db; do
	echo $f
	gemini bcolz_index $f;
done
