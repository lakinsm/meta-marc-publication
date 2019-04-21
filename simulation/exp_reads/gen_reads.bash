#!/bin/bash
cd ..
basePath=`pwd`
cd exp_reads

db="$basePath/megares_database_v1.01.fasta"
fr="$basePath/bin/fake_gene"

Genes=( 0 3 5 )
Y=( 0 0.1 0.25 0.5 0.7 )
Z=( 5 25 50 )
Contiguous=( "true" "false" )

for g in "${Genes[@]}"
do
    for y in "${Y[@]}"
    do
	    for z in "${Z[@]}"
	    do
	        for c in "${C[@]}"
                do
	            echo "Generating the synthetic gene..."
	            cmd="$fr $db $y $z $g $c"
	            echo $cmd
	            $cmd
	            outFile="fake-gene=$g-$z-$y-$c.fa"
	            echo "Generating fake reads..."
	            cov='10'
	            nReads=`wc -m $outFile | awk '{print $1}'`
	            nReads=`echo "$nReads / 100 * $cov" | bc`
	            echo "nReads=$nReads"
	            readFileBase="reads-gene=$g-$z-$y-$c"
	            cmd="$basePath/ecoli/BEAR/scripts/generate_reads.py -r $outFile -a abund.txt -o $readFileBase -t $nReads -l 100 -i 100 -s 1 -d"
	            echo $cmd
	            $cmd
	            readFileOne="${readFileBase}.1.fasta"
	            readFileTwo="${readFileBase}.2.fasta"
	            echo "Concatenating reads with ecoli..."
	            cat $basePath/ecoli/ecoli.1.fasta $readFileOne > reads-$g-$z-$y-$c.1.fna
	            cat $basePath/ecoli/ecoli.2.fasta $readFileTwo > reads-$g-$z-$y-$c.2.fna
	        done
	    done
    done
done
rm *.fa
rm *.fasta

