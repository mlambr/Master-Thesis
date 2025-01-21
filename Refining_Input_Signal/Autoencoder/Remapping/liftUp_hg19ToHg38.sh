#!/bin/bash
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4
#SBATCH --array=2-15
#SBATCH --output=liftUp_hg19ToHg38_%A-%a.out    # Standard output and error log


## job name
#SBATCH --job-name="liftUp_hg19ToHg38"

## call before running this script
# conda activate cfdna


BAM_LIST=($(<dsDNA_list.txt))
BAM=${BAM_LIST[${SLURM_ARRAY_TASK_ID}]}

echo "CrossMap.py bam -a /cluster/work/medinfmk/ref_genomes/hg19ToHg38.over.chain.gz --chromid l $BAM hg38/${BAM/.bam/}"
CrossMap.py bam -a /cluster/work/medinfmk/ref_genomes/hg19ToHg38.over.chain.gz --chromid l $BAM hg38/${BAM/.bam/}

#MAXTHREADS=5
#THREADS=$MAXTHREADS
#PID=
#while read bam; do
#for bam in *.bam; do
#	if [[ THREADS == 0 && ]]; then
#		echo "wait $PID
#		wait $PID
#		THREADS=$MAXTHREADS
#		fi
#	echo "CrossMap.py bam -a /cluster/work/medinfmk/ref_genomes/hg19ToHg38.over.chain.gz --chromid l $bam hg38/${bam/.bam/}" & 
#	# CrossMap.py bam -a /cluster/work/medinfmk/ref_genomes/hg19ToHg38.over.chain.gz --chromid l $bam hg38/${bam/.bam/}
#	PID=$!
#	let THREADS--
#done < dsDNA_list.txt
#done
