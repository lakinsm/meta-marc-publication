#!/usr/bin/env bash

if [ -n "$(ls -A ${PWD}/simulated_data)" ]; then
    rm "simulated_data/*"
fi

for insertion in 0.1 0.25 0.5 0.7; do
    for copy_num in 5 25 50; do
        for contiguous in "true" "false"; do
        scripts/insert_sequences.py megares_database_v1.01.fasta 300 ${copy_num} ${insertion} ${contiguous} simulated_data
        done
    done
done

