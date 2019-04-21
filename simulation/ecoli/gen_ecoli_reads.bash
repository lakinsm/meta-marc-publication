#!/bin/bash
cmd="cat *.fna > ecoli.fna"
echo $cmd
cat *.fna > ecoli.fna

ecoli='4600000'
cov='10'
nReads=`echo "$ecoli / 100 * $cov" | bc `
echo "nReads = $nReads"
cmd="./BEAR/scripts/generate_reads.py -r ./ecoli.fna -a abund.txt -o ecoli -t $nReads -l 100 -i 100 -s 1 -d"
echo $cmd
$cmd

rm ecoli.fna