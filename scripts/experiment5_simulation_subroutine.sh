#!/usr/bin/env bash

parameter_string=$( basename $1 | sed -r 's/\.fasta$//' )

echo "$parameter_string"

mkdir temp/${parameter_string}
mkdir output/${parameter_string}

# Data generation
art_illumina -ss HS25 -i $1 -l 150 -f 10 -m 250 -s 10 -o "temp/${parameter_string}/test_seqs_simulated" -rs 154 -p -na -q
cat ecoli_data/all_ecoli_simulated1.fq >> temp/${parameter_string}/test_seqs_simulated1.fq
cat ecoli_data/all_ecoli_simulated2.fq >> temp/${parameter_string}/test_seqs_simulated2.fq

# Meta-MARC HTS Reads
grep -A 1 "^@" temp/${parameter_string}/test_seqs_simulated1.fq | grep -v "\-\-" | sed 's/^@/>/' > temp/${parameter_string}/test_seqs_simulated.fasta
grep -A 1 "^@" temp/${parameter_string}/test_seqs_simulated2.fq | grep -v "\-\-" | sed 's/^@/>/' >> temp/${parameter_string}/test_seqs_simulated.fasta

nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_hts_outfileI.tblout.scan mmarc_hmms/mmarc_groupI.hmm temp/${parameter_string}/test_seqs_simulated.fasta > /dev/null
nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_hts_outfileII.tblout.scan mmarc_hmms/mmarc_groupII.hmm temp/${parameter_string}/test_seqs_simulated.fasta > /dev/null
nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_hts_outfileIII.tblout.scan mmarc_hmms/mmarc_groupIII.hmm temp/${parameter_string}/test_seqs_simulated.fasta > /dev/null

# Alignment
java -jar /home/lakinsm/bin/trimmomatic-0.36.jar PE -threads 2 -phred33 temp/${parameter_string}/test_seqs_simulated1.fq temp/${parameter_string}/test_seqs_simulated2.fq temp/${parameter_string}/trimmed_1P.fastq temp/${parameter_string}/trimmed_1U.fastq temp/${parameter_string}/trimmed_2P.fastq temp/${parameter_string}/trimmed_2U.fastq ILLUMINACLIP:/home/lakinsm/amr_tools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
bwa mem -t 2 databases/mmarc_groupI_seqs.fasta temp/${parameter_string}/trimmed_1P.fastq temp/${parameter_string}/trimmed_2P.fastq > output/${parameter_string}/alignment_I.sam
bwa mem -t 2 databases/mmarc_groupII_seqs.fasta temp/${parameter_string}/trimmed_1P.fastq temp/${parameter_string}/trimmed_2P.fastq > output/${parameter_string}/alignment_II.sam
bwa mem -t 2 databases/mmarc_groupIII_seqs.fasta temp/${parameter_string}/trimmed_1P.fastq temp/${parameter_string}/trimmed_2P.fastq > output/${parameter_string}/alignment_III.sam

# Base assembly
fq2fa --merge --filter temp/${parameter_string}/test_seqs_simulated1.fq temp/${parameter_string}/test_seqs_simulated2.fq temp/${parameter_string}/interleavened.fasta
max_len=$( grep -v "^>" temp/${parameter_string}/interleavened.fasta | wc -L )

if (( "$max_len" <= "128" )); then
	idba_ud --num_threads 2 -r temp/${parameter_string}/interleavened.fasta -o temp/${parameter_string}/idba
else
	idba_ud --num_threads 2 -l temp/${parameter_string}/interleavened.fasta -o temp/${parameter_string}/idba
fi

cp temp/${parameter_string}/idba/contig.fa contigs/${parameter_string}_contig.fa
mv temp/${parameter_string}/idba/contig.fa temp/${parameter_string}/contig.fa

rm temp/${parameter_string}/idba/*
rmdir temp/${parameter_string}/idba

# Alignment to contigs
bwa index temp/${parameter_string}/contig.fa
bwa mem -t 2 temp/${parameter_string}/contig.fa temp/${parameter_string}/trimmed_1P.fastq temp/${parameter_string}/trimmed_2P.fastq > output/${parameter_string}/contig_alignments.sam

# Resfams contig annotation
gmhmmp -m /home/lakinsm/MetaGeneMark_linux_64/MetaGeneMark/MetaGeneMark_v1.mod -A temp/${parameter_string}/protein.fasta temp/${parameter_string}/contig.fa
cat temp/${parameter_string}/protein.fasta | grep -v -e "^$" > temp/${parameter_string}/protein_clean.fasta
hmmsearch --notextw --cpu 2 --tblout output/${parameter_string}/resfams_contig.tblout.scan resfams_hmms/Resfams-full.hmm temp/${parameter_string}/protein_clean.fasta > /dev/null

# Meta-MARC contig annotation
nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_contig_I.tblout.scan mmarc_hmms/mmarc_groupI.hmm temp/${parameter_string}/contig.fa > /dev/null
nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_contig_II.tblout.scan mmarc_hmms/mmarc_groupII.hmm temp/${parameter_string}/contig.fa > /dev/null
nhmmer --dna --notextw --cpu 2 --tblout output/${parameter_string}/mmarc_contig_III.tblout.scan mmarc_hmms/mmarc_groupIII.hmm temp/${parameter_string}/contig.fa > /dev/null

rm temp/${parameter_string}/*
rmdir temp/${parameter_string}

