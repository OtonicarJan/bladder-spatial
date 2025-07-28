#!/bin/bash


# LSF directives
#BSUB -J "ichorCNA[1-68]"
#BSUB -e logs/cna-%J-%I.err
#BSUB -o logs/cna-%J-%I.out
#BSUB -q medium
#BSUB -R "rusage[mem=20G]"
#BSUB -P ichorCNA

module load R/4.3.0


FASTQ_PATHS="chromopthripsis/bladder_cancer/liquid_biopsy/fastq_paths.txt"

mkdir -p chromopthripsis/bladder_cancer/liquid_biopsy/bam_files
cd chromopthripsis/bladder_cancer/liquid_biopsy/bam_files

# Extract sample name and fastq files and run a job
sample_num=${LSB_JOBINDEX}
sample_id=$(awk -v sample_num="$sample_num" 'NR==sample_num {print $1}' $FASTQ_PATHS)

cd $sample_id

echo "Working on sample $sample_id"

Rscript chromopthripsis/j462r/tools/ichorCNA/scripts/runIchorCNA.R --id 'ichorCNA_'$sample_id \
--WIG 'chromopthripsis/bladder_cancer/liquid_biopsy/bam_files/'$sample_id'/'$sample_id'.size_exclusion.500.kb.wig' --ploidy "c(2,3,4)" --normal "c(0.7,0.8,0.9,0.95,0.99,0.995,0.999,0.9999)" --maxCN 6 \
                                --gcWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.gc.seg \
                                --mapWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.map.seg \
                                --centromere chromopthripsis/j462r/tools/ichorCNA/inst/extdata/GRCh38.GCA_000001405.26_centromere_acen.txt \
                                --normalPanel chromopthripsis/bladder_cancer/liquid_biopsy/panel_of_normal_median.rds \
                                --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --plotYLim "c(-1,1)" \
                                --estimateNormal True --estimatePloidy True --estimateScPrevalence False \
                                --scStates "c()" --txnE 0.9999 --txnStrength 10000 --genomeBuild "hg38" --genomeStyle "UCSC"

Rscript chromopthripsis/j462r/tools/ichorCNA/scripts/runIchorCNA.R --id 'ichorCNA_all_'$sample_id \
--WIG 'chromopthripsis/bladder_cancer/liquid_biopsy/bam_files/'$sample_id'/'$sample_id'.500.kb.wig' --ploidy "c(2,3,4)" --normal "c(0.7,0.8,0.9,0.95,0.99,0.995,0.999,0.9999)" --maxCN 6 \
                                --gcWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.gc.seg \
                                --mapWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.map.seg \
                                --centromere chromopthripsis/j462r/tools/ichorCNA/inst/extdata/GRCh38.GCA_000001405.26_centromere_acen.txt \
                                --normalPanel chromopthripsis/bladder_cancer/liquid_biopsy/panel_of_normal_all_fragments_median.rds \
                                --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
                                --estimateNormal True --estimatePloidy True --estimateScPrevalence False --plotYLim "c(-1,1)" \
                                --scStates "c()" --txnE 0.9999 --txnStrength 10000 --genomeBuild "hg38" --genomeStyle "UCSC"