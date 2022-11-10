#!/bin/bash

REF=./ref/
GENCODE=./gencode/
FASTA='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz'
GTF='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz'

if [ -d "$REF" ] || [ -d "$GENCODE" ];
then
    echo "$REF directory exists."
    echo "$GENCODE directory exists"
else
	echo "$REF and $GENCODE directory does not exist...setting up references"
	mkdir $REF $GENCODE
	wget $FASTA -P $GENCODE; wget $GTF -P $GENCODE
	gunzip -f ./gencode/GRCh38.primary_assembly.genome.fa.gz
	gunzip -f ./gencode/gencode.v42.primary_assembly.annotation.gtf.gz
fi
