#!/usr/bin/env bash

num_simulated=$( scripts/count_num_simulated_from_sam.py ${1}/alignment_I.sam )
params=$( basename $1 )

# Accumulate alignment results
scripts/count_num_aligned_from_sam.py ${1}/alignment_I.sam ${1}/alignment_II.sam ${1}/alignment_III.sam megares_annotations_v1.01.csv $num_simulated > results/${params}.csv

# Accumulate Meta-MARC HTS results
scripts/count_num_classified_mmarc_hts.py ${1}/mmarc_hts_outfileI.tblout.scan ${1}/mmarc_hts_outfileII.tblout.scan ${1}/mmarc_hts_outfileIII.tblout.scan megares_annotations_v1.01.csv mmarc_model_annotations.tsv mmarc_model_members.csv $num_simulated >> results/${params}.csv

# Accumulate Meta-MARC Assembly results
cat ${1}/contig_alignments.sam | scripts/count_num_classified_mmarc_assembly.py ${1}/mmarc_contig_I.tblout.scan ${1}/mmarc_contig_II.tblout.scan ${1}/mmarc_contig_III.tblout.scan megares_annotations_v1.01.csv mmarc_model_annotations.tsv mmarc_model_members.csv $num_simulated >> results/${params}.csv

# Accumulate Resfams results
cat ${1}/contig_alignments.sam | scripts/count_num_classified_resfams.py ${1}/resfams_contig.tblout.scan megares_annotations_v1.01.csv resfams_metadata.txt $num_simulated  >> results/${params}.csv

