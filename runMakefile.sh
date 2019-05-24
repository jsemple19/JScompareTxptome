#!/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
#SBATCH --time=0-24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=4G
#SBATCH --array=6-96
#SBATCH --output=/home/ubelix/izb/semple/output/slurm-%j.out
#SBATCH --error=/home/ubelix/izb/semple/error/slurm-%j.err

module add vital-it
module add UHTS/Quality_control/fastqc/0.11.5      #fastqc
module add UHTS/Quality_control/cutadapt/1.13     #cutadapt
module add UHTS/Analysis/trimmomatic/0.36;
module add UHTS/Aligner/STAR/2.6.0c;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/salmon/0.11.2;


fastqFileList=(`cat ../fastqList.txt`)


let i=${SLURM_ARRAY_TASK_ID}-1
#fastqFile=${fastqFileList[$i]}

make all fastqFile=${fastqFileList[$i]}


