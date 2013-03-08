#!/bin/bash
usage() { 
    echo "usage: sh $0 <input vcf> <num. processors> <gemini db name>"  
}

if [ -z $3 ]
then
	usage
	exit
fi

VCF=$1
PROCS=$2
DB=$3

# BGZIP and grabix the VCF so that we can grabix arbitrary chunks.
echo "bgzipping the VCF"
bgzip $VCF

echo "making a grabix index of the bgzipped VCF for chunking"
grabix index $VCF".gz"

# how many lines are in the file?
NUM_LINES=$(awk 'NR==2' $VCF".gz.gbi")
CHUNKSIZE=$(($NUM_LINES / $PROCS))
GRABIXFILE=$VCF".gz"
CHUNKNUM=1
CHUNKFILES=()

echo "chunking VCF and loading each chunk"

MAX=$(($NUM_LINES-1))
for i in `seq 0 $CHUNKSIZE $MAX`
do
    from=$(($i+1))
    to=$(($i+$CHUNKSIZE))
    
    if [ $to -gt $NUM_LINES ]
    then
      to=$NUM_LINES
    fi
    
    grabix grab $GRABIXFILE $from $to > $VCF".chunk"$CHUNKNUM
    gemini load_chunk -v $VCF".chunk"$CHUNKNUM -t snpEff $VCF".chunk"$CHUNKNUM".db" -o $from &
    
    CHUNKFILES+=($VCF".chunk"$CHUNKNUM".db")
    CHUNKNUM=$(($CHUNKNUM+1))
    
done

wait

echo "merging chunks"

CHUNKLIST=""
for chunkfile in "${CHUNKFILES[@]}"
do
    CHUNK=" --chunkdb "$chunkfile
    CHUNKLIST=$CHUNKLIST$CHUNK
done

gemini merge_chunks $CHUNKLIST --db foo.db

# cleanup temp files
rm $VCF".chunk"*

echo "complete"

