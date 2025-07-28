#!/bin/bash

# an example of spaceranger run
# run it with 'bsub < run_spaceranger.sh'

#BSUB -J B4
#BSUB -R "rusage[mem=150G]"
#BSUB -n 16
#BSUB -q long

module load spaceranger/2.0.1

SAMPLE='B4'
cd spatial_transcriptomics/spaceranger_results


spaceranger count --id=${SAMPLE} \
--transcriptome=../reference/refdata-gex-GRCh38-2020-A \
--fastqs=../raw_fastqs/run230509_A01382_0373_BH5FHWDRX3/${SAMPLE} \
--image=../images/Spatial_experiment_1_images/export/B420_FFPE_run1_20230418.vsi.Collection/B420_FFPE_run1_20230418_20x_BF_EFI_02.tif \
--probe-set=/software/spaceranger/2.0.1/probe_sets/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \
--slide=V12U06-071 \
--area=B1 \
--localcores=16 \
--localmem=150
