#!/bin/bash

# Workflow to process liquid biopsy samples
# Additionally, BAM files are subset to only reads with insert size between 90 and 150bp
# as done by Sottoriva lab https://www.nature.com/articles/s43018-024-00787-0#Sec11 

# LSF directives
#BSUB -J "liquid_biopsy[1-68]"
#BSUB -e logs/T-%J-%I.err
#BSUB -o logs/T-%J-%I.out
#BSUB -q long
#BSUB -n 8
#BSUB -R "rusage[mem=50G]"
#BSUB -P liquid_biopsy

module load samtools
module load picard/2.25.1
module load bwa
module load gcc/11.1.0
module load qualimap

export JAVA_OPTS="-Xmx40G"


FASTQ_PATHS="chromopthripsis/bladder_cancer/liquid_biopsy/fastq_paths.txt"
REFERENCE="chromopthripsis/j462r/references/hg38_ONT/1KG_ONT_VIENNA_hg38.fa"

mkdir -p chromopthripsis/bladder_cancer/liquid_biopsy/bam_files
cd chromopthripsis/bladder_cancer/liquid_biopsy/bam_files

# Extract sample name and fastq files and run a job
sample_num=${LSB_JOBINDEX}
sample_id=$(awk -v sample_num="$sample_num" 'NR==sample_num {print $1}' $FASTQ_PATHS)
forward_read=$(awk -v sample_num="$sample_num" 'NR==sample_num {print $2}' $FASTQ_PATHS)
reverse_read=$(awk -v sample_num="$sample_num" 'NR==sample_num {print $3}' $FASTQ_PATHS)

echo 'Sample number = '$sample_num
echo 'Sample ID = '$sample_id

mkdir -p $sample_id
cd $sample_id

ALIGNED_BAM=$sample_id'.aligned.bam'
CELL_BAM=$sample_id'.sorted.bam'
DUPLICATES=$sample_id'.sorted.mdup.bam'
CLEAN_BAM=$sample_id'.filtered.bam'
CLEAN_SORT_BAM=$sample_id'.filtered.sorted.bam'
MDUP_METRICS=$sample_id'.mdup.metrics.txt'
BAM_FILE_90_150=$sample_id'.size_exclusion.bam'
BAM_FILE_90_150_SORTED=$sample_id'.sorted.size_exclusion.bam'
BAM_FILE_151_220=$sample_id'.size_exclusion_large.bam'
BAM_FILE_151_220_SORTED=$sample_id'.sorted.size_exclusion_large.bam'

if [ -f "$CLEAN_SORT_BAM" ]; then
    echo "$CLEAN_SORT_BAM exists."
else 
    echo "$CLEAN_SORT_BAM does not exist. Will be created"
    
    if [ -f "$CELL_BAM" ]; then
        echo "$CELL_BAM exists."
    else
        echo "Running alignment for $sample_id"
        bwa mem -M -t 8 $REFERENCE $forward_read $reverse_read | samtools view -bS -o $ALIGNED_BAM
        samtools sort -@ 8 -o $CELL_BAM $ALIGNED_BAM
        samtools index -@ 8 $CELL_BAM
    fi
    
    if [ -f "$DUPLICATES" ]; then
        echo "$DUPLICATES exists."
    else
        echo "Marking duplicates and filtering BAM file"
        picard.sh MarkDuplicates I=$CELL_BAM O=$DUPLICATES M=$MDUP_METRICS
        
    fi
    
    samtools view -h -b -F 3844 -q 30 $DUPLICATES > $CLEAN_BAM
    samtools sort -@ 8 $CLEAN_BAM -o $CLEAN_SORT_BAM
    samtools index -@ 8 $CLEAN_SORT_BAM
    samtools flagstat $DUPLICATES > 'flagstat.'$sample_id'.original.txt'
    samtools flagstat $CLEAN_SORT_BAM > 'flagstat.'$sample_id'.filtered.txt'

    
    rm $ALIGNED_BAM
    rm $DUPLICATES 
    rm $CLEAN_BAM
    rm $MDUP_METRICS
fi


# Extracting only the pairs of reads with insert size between 90 and 150bp
if [ -f "$BAM_FILE_90_150_SORTED" ]; then
    echo "$BAM_FILE_90_150_SORTED exists."
else
    echo "Subsetting cleaned BAM file to insert sizes between 90 and 150bp"
    samtools view -h $CLEAN_SORT_BAM | \
      awk 'substr($0,1,1)=="@" || ($9>=90 && $9<=150) || ($9<=-90 && $9>=-150)' | \
      samtools view -b > $BAM_FILE_90_150
    samtools sort -@ 8 $BAM_FILE_90_150 -o $BAM_FILE_90_150_SORTED
    samtools index -@ 8 $BAM_FILE_90_150_SORTED
    samtools flagstat $BAM_FILE_90_150_SORTED > 'flagstat.'$sample_id'.size_excl.txt'
    
    rm $BAM_FILE_90_150
fi

# Extracting only the pairs of reads with insert size between 151 and 220bp
if [ -f "$BAM_FILE_151_220_SORTED" ]; then
    echo "$BAM_FILE_151_220_SORTED exists."
else
    echo "Subsetting cleaned BAM file to insert sizes between 90 and 150bp"
    samtools view -h $CLEAN_SORT_BAM | \
      awk 'substr($0,1,1)=="@" || ($9>=151 && $9<=220) || ($9<=-151 && $9>=-220)' | \
      samtools view -b > $BAM_FILE_151_220
    samtools sort -@ 8 $BAM_FILE_151_220 -o $BAM_FILE_151_220_SORTED
    samtools index -@ 8 $BAM_FILE_151_220_SORTED
    samtools flagstat $BAM_FILE_151_220_SORTED > 'flagstat.'$sample_id'.size_excl_large.txt'
    
    rm $BAM_FILE_151_220
fi


# Run readCounter for different bin sizes
bins_sizes=(100 250 500 1000)
for bin_kb in "${bins_sizes[@]}"; do
    bin=${bin_kb}000
    SEGMENT_FILE="$sample_id.$bin_kb.kb.wig"
    SEGMENT_FILE_SIZE="$sample_id.size_exclusion.$bin_kb.kb.wig"
    SEGMENT_FILE_SIZE_LARGE="$sample_id.size_exclusion_large.$bin_kb.kb.wig"
    if [ -f "$SEGMENT_FILE" ] && [ -s "$SEGMENT_FILE" ]; then
        echo "$SEGMENT_FILE exists and is not empty."
    else 
        echo "$SEGMENT_FILE does not exist. Will be created"
        chromopthripsis/j462r/tools/hmmcopy_utils/bin/readCounter -q 30 -w ${bin} -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $CLEAN_SORT_BAM > $SEGMENT_FILE
    fi
    
    if [ -f "$SEGMENT_FILE_SIZE" ] && [ -s "$SEGMENT_FILE_SIZE" ]; then
        echo "$SEGMENT_FILE_SIZE exists and is not empty."
    else 
        echo "$SEGMENT_FILE_SIZE does not exist. Will be created"
        chromopthripsis/j462r/tools/hmmcopy_utils/bin/readCounter -q 30 -w ${bin} -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $BAM_FILE_90_150_SORTED > $SEGMENT_FILE_SIZE
    fi
    
    if [ -f "$SEGMENT_FILE_SIZE_LARGE" ] && [ -s "$SEGMENT_FILE_SIZE_LARGE" ]; then
        echo "$SEGMENT_FILE_SIZE_LARGE exists and is not empty."
    else 
        echo "$SEGMENT_FILE_SIZE_LARGE does not exist. Will be created"
        chromopthripsis/j462r/tools/hmmcopy_utils/bin/readCounter -q 30 -w ${bin} -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $BAM_FILE_151_220_SORTED > $SEGMENT_FILE_SIZE_LARGE
    fi
    
done

if [ -d "stats" ]; then
    echo "Stats for sample $sample_id exist"
else
    qualimap bamqc -bam $CLEAN_SORT_BAM --java-mem-size=38G -nt 6 -outdir stats
fi