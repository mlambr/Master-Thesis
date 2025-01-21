#!/bin/bash

#SBATCH --job-name=M1_7
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --output=/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/run/training_model.out


{
    start_time=$(date +%s)

    # Show memory usage with sinfo
    python3 -u /cluster/work/medinfmk/cfDNA-Snyder/results/Maria/run/Model_corrected.py 

    end_time=$(date +%s)
    execution_time=$((end_time - start_time))

    echo "Execution time: $execution_time seconds"
} >> /cluster/work/medinfmk/cfDNA-Snyder/results/Maria/run/training_model.out 2>&1
