
db=extended_ped.db
gemini bcolz_index $db
t=1

for filter in \
	'(gt_types.M10475 == HET) and (gt_types.M10478 == HET)' \
	'(gt_types.M10478 == HET or gt_types.M10475 == HET)' \
	'(gt_types).(phenotype == 1).(==HOM_REF).(all) and ((gt_types).(phenotype==2).(==HET).(any) or (gt_types).(phenotype==2).(==HOM_ALT).(any))' \
	'(gt_types).(phenotype == 1).(==HOM_REF).(all) and gt_types.M10475 == HET' \
	'(gt_types).(phenotype == 1).(==HOM_REF).(all) and gt_types.M10475 == HET' \
	'(gt_types).(phenotype==2).(!=HOM_REF).(all) and (gt_depths).(phenotype==2).(<10).(none)' \
	'(gt_types).(*).(!=HOM_REF).(all) and (gt_depths).(phenotype==2).(<10).(any)' \
	'(gt_types).(*).(!=HOM_REF).(any) or (gt_depths).(*).(==HOM_REF).(any)' \
	'(gt_depths).(phenotype==2).(<10).(count <= 2)' \
	; do
	# NOTE that count queries must be on their own.
	echo "bcolz-test.t$t"
	t=$((t + 1))
	gemini query -q "select chrom, start, end from variants" --gt-filter "$filter" $db > exp
	gemini query -q "select chrom, start, end from variants" --gt-filter "$filter" $db --use-bcolz > obs
	check obs exp
	echo ""
done
