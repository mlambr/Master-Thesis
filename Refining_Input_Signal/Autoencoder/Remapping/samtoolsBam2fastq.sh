#!/bin/bash
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=8
##SBATCH --array=4 # if you want a specific sample only
#SBATCH --output=samtoolsBam2fastq_%A-%a.out    # Standard output and error log


## job name
#SBATCH --job-name="samtoolsBam2fastq"

## call before running this script
# conda activate cfdna

# enable_modules
module load samtools

# mkdir fastq_st
	
BAM_LIST=($(<dsDNA_list.txt))
BAM=${BAM_LIST[${SLURM_ARRAY_TASK_ID}]}

SAMPLE=${BAM/.bam/}
echo "samtools fastq -@ 4 $BAM -1 fastq_st/${SAMPLE}_R1.fastq.gz -2 fastq_st/${SAMPLE}_R2.fastq.gz -0 /dev/null -s /dev/null -n"

if [ ! -f ${SAMPLE}.sorted.bam ]
then
	echo "samtools sort -@ 4 -n $BAM -o ${SAMPLE}.sorted.bam"
	samtools sort -@ 4 -n $BAM -o ${SAMPLE}.sorted.bam
fi

samtools fastq -@ 4 ${SAMPLE}.sorted.bam -1 fastq_st/${SAMPLE}_R1.fastq.gz -2 fastq_st/${SAMPLE}_R2.fastq.gz -0 /dev/null -s /dev/null -n
