#!/usr/bin/bash

FASTQ_FILE=$1

if [[ $FASTQ_FILE == *.gz ]]; then
    tar -zxvf $FASTQ_FILE
    if [[ $(ls *fastq.gz) == *.gz ]]; then
        gunzip -d *fastq.gz
    fi
    sleep 30
else
    tar -xvf $FASTQ_FILE; gunzip -d *fastq.gz
    sleep 30
fi

