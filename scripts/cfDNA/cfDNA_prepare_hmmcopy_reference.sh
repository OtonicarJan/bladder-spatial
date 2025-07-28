#!/bin/sh
# LSF directivesy
#BSUB -J "mappability[1-4]"
#BSUB -e logs/T-%J-mappability.err
#BSUB -o logs/T-%J-mappability.out
#BSUB -q verylong
#BSUB -n 5
#BSUB -R "rusage[mem=80GB]"
#BSUB -P mappability

echo + `date` $LSB_JOBNAME started on $HOSTNAME jobID=$LSB_JOBID 

module load gcc/11.1.0
module load bowtie/1.3.0 

bins=(100 250 500 1000)

bin_kb=${bins[$((LSB_JOBINDEX-1))]}
#bin_kb=1000
bin=${bin_kb}000
REF='chromopthripsis/j462r/references/hg38_ONT/1KG_ONT_VIENNA_hg38.fa'

cd chromopthripsis/j462r/references/hg38_ONT/mappability_gc

## Mappability
if [ -f "hg38.150len.fa.map.bw" ]; then
    echo "hg38.150len.fa.map.bw  exists."
else 
    echo "hg38.150len.fa.map.bw  does not exist. Will be created"
    chromopthripsis/j462r/tools/hmmcopy_utils/util/mappability/generateMap.pl --window 150 -o hg38.150len.fa.map.bw  -b $REF
    chromopthripsis/j462r/tools/hmmcopy_utils/util/mappability/generateMap.pl --window 150 -o hg38.150len.fa.map.bw  $REF
fi


chromopthripsis/j462r/tools/hmmcopy_utils/bin/mapCounter -w ${bin}  -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY chromopthripsis/j462r/references/hg38_ONT/mappability_gc/hg38.150len.fa.map.bw > hg38.${bin_kb}kb.map.seg

## GC bias
chromopthripsis/j462r/tools/hmmcopy_utils/bin/gcCounter -w ${bin}  -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $REF > hg38.${bin_kb}kb.gc.seg
