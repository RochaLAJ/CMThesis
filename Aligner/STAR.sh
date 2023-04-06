#!/bin/bash

STAR --runMode genomeGenerate \
--genomeDir ./ref/ \
--genomeFastaFiles ./gencode/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ./gencode/gencode.v42.primary_assembly.annotation.gtf \
--genomeSAsparseD 2 \
--genomeSAindexNbases 2 \ 
--runThreadN 16 \
--genomeChrBinNbits 6  \
--limitGenomeGenerateRAM 20000000000


STAR --runMode alignReads \
--outSAMtype BAM Unsorted SortedByCoordinate \
--runThreadN 16 \
--genomeDir ./ref/ \
--outFilterMultimapNmax 50 \
--readFilesIn read1.fastq \
--quantMode GeneCounts TranscriptomeSAM \
--outReadsUnmapped Fastx
