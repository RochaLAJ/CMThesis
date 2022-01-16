#!/bin/bash

STAR --runMode genomeGenerate --genomeDir Genome --ge gnomeFastaFiles /star_index/ref.fa  \ 
         --genomeSAsparseD 2 --genomeSAindexNbases 2 --runThreadN 16 --genomeChrBinNbits 6

STAR --runMode alignReads --outSAMtype BAM Unsorted SortedByCoordinte --runThreadN 16 --genomeDir Genome --outFilterMultimapNmax 50 
       \ --readFilesIn /chip/bam49.fasta --quantMode GeneCounts
       \  TranscriptomeSAM --outTmpKeep /run/media/lemt/Backup/bam49 --outReadsUnmapped Fastx
