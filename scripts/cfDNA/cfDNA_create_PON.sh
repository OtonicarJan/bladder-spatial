#!/bin/bash

# LSF directives
#BSUB -q medium
#BSUB -R "rusage[mem=20G]"
#BSUB -P ichorCNA

module load R/4.3.0

cd chromothripsis/bladder_cancer/liquid_biopsy/panel_of_normals

# input is just .wig file for the benign sample
Rscript chromopthripsis/j462r/tools/ichorCNA/scripts/createPanelOfNormals.R \
--filelist chromopthripsis/bladder_cancer/liquid_biopsy/panel_of_normals/PON_wig_500kb.txt \
--gcWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.gc.seg \
--mapWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.map.seg \
--centromere chromopthripsis/j462r/tools/ichorCNA/inst/extdata/GRCh38.GCA_000001405.26_centromere_acen.txt \
--outfile panel_of_normal_all_fragments_median

Rscript chromopthripsis/j462r/tools/ichorCNA/scripts/createPanelOfNormals.R \
--filelist chromopthripsis/bladder_cancer/liquid_biopsy/panel_of_normals/PON_wig_short_500kb.txt \
--gcWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.gc.seg \
--mapWig chromothripsis/j462r/references/hg38_ONT/mappability_gc/hg38.500kb.map.seg \
--centromere chromopthripsis/j462r/tools/ichorCNA/inst/extdata/GRCh38.GCA_000001405.26_centromere_acen.txt \
--outfile panel_of_normal_median