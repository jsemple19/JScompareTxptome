#! /usr/bin/bash
#SBATCH --mail-user=jennifer.semple@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="RNAseq"
#SBATCH --time=0-24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1

module add vital-it
#module add UHTS/Quality_control/fastqc/0.11.5
module add UHTS/Quality_control/fastqc/0.11.7;
#module add UHTS/Quality_control/cutadapt/1.13
module add UHTS/Quality_control/cutadapt/2.5;
#module add UHTS/Analysis/trimmomatic/0.36;
#module add UHTS/Aligner/STAR/2.6.0c;
module add UHTS/Aligner/STAR/2.7.3a;
module add UHTS/Analysis/samtools/1.10;
#module add UHTS/Analysis/samtools/1.8;
#module add UHTS/Analysis/salmon/0.11.2;
export SALMON_SING="singularity exec /software/singularity/containers/salmon-1.2.1-1.ubuntu18.sif"

fastqFileList=./fastqList.txt
fastqFile=(`cut -f1 $fastqFileList`)
sampleName=(`cut -f2 $fastqFileList`)
repeatNum=(`cut -f3 $fastqFileList`)

i=${SLURM_ARRAY_TASK_ID}

#####################
# mapping RNAseq data
#####################

###############################
########### VARIABLES #########
###############################


fastqFile=${fastqFile[$i]}
sampleName=${sampleName[$i]}
repeatNum=${repeatNum[$i]}
nThreads=${SLURM_CPUS_PER_TASK}

genomeVer=WS260
genomeFile=${HOME}/genomeVer/${genomeVer}/sequence/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
GENOME_DIR=`dirname ${genomeFile}`
# annotFile=/home/ubelix/izb/semple/genomeVer/${genomeVer}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gff3.gz
# gunzip $annotFile
# annotFile=${annotFile%.gz}
# use gffread from cufflinks to convert gff to gtf
# module add UHTS/Assembler/cufflinks/2.2.1
# b=(`basename -s .gff3 ${annotFile}`)
# gffread $annotFile -T -o ${annotFile%/*}/${b}.gtf
# need to remove wierd exons:
# grep WormBase c_elegans.PRJNA13758.${genomeVer}.annotations.gtf > c_elegans.PRJNA13758.WS260.annotations1.gtf
# mv c_elegans.PRJNA13758.${genomeVer}.annotations1.gtf c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
annotFile=${HOME}/genomeVer/${genomeVer}/annotation/c_elegans.PRJNA13758.${genomeVer}.annotations.gtf
#mRNAseqFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.mRNA_transcripts.fa.gz
mRNAindex=/home/ubelix/izb/semple/genomeVer/${genomeVer}/sequence/${genomeVer}_mRNA_index
ncRNAindex=/home/ubelix/izb/semple/genomeVer/${genomeVer}/sequence/${genomeVer}_ncRNA_index
pseudoIndex=/home/ubelix/izb/semple/genomeVer/${genomeVer}/sequence/${genomeVer}_pseudogenic_index
tnIndex=/home/ubelix/izb/semple/genomeVer/${genomeVer}/sequence/${genomeVer}_transposon_index


WORK_DIR=$PWD
QC_DIR=${WORK_DIR}/qc

#FASTQ_DIR=`dirname ${fastqFile}`
baseName=${sampleName}_${repeatNum}

########################################################
### get initial read stats                            ##
########################################################
#
##run fastqc on sequences
#mkdir -p ${WORK_DIR}/qc/rawData
#fastqc ${fastqFile} -o ${WORK_DIR}/qc/rawData 
#	
########################################################
### trim adaptors with cutadapt                       ##
########################################################
#
## use cutadapt to trim
#mkdir -p cutadapt
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o cutadapt/${baseName}.fastq.gz  ${fastqFile}
#
##redo fastQC on trimmed reads
#mkdir -p ${WORK_DIR}/qc/cutadapt
#fastqc cutadapt/${baseName}.fastq.gz -o ${WORK_DIR}/qc/cutadapt

		
#######################################################
## Align to genome with STAR                         ##
#######################################################


# index genome
#echo "indexing genome..."
#STAR --runMode genomeGenerate --genomeDir ${genomeFile%/*} --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotFile} --runThreadN $nThreads

# align to genome
echo "aligning to genome..."
mkdir -p ${WORK_DIR}/bam
STAR --genomeDir ${GENOME_DIR}  --readFilesIn ${WORK_DIR}/cutadapt/${baseName}.fastq.gz --readFilesCommand zcat --outFileNamePrefix ${WORK_DIR}/bam/${baseName} --runThreadN $nThreads --alignIntronMax 500 --quantMode GeneCounts

# convert sam to bam
samtools view -b ${WORK_DIR}/bam/${baseName}.sam -o ${WORK_DIR}/bam/${baseName}.bam
rm ${WORK_DIR}/bam/${baseName}.sam

#######################################################
## Count reads with Salmon                           ##
#######################################################

# prepare reference genome
# See separate script indexTxpts4salmon.sh

######## NOTE: I do not have estimates for --fldMean and --fldSD as i have no access to the bioanalyser files #########
#######  therefore the quantification based on the effective transcript length will be wrong!!!! ######

# quantify mRNA transcripts
${SALMON_SING} salmon quant -i ${mRNAindex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz -o ${WORK_DIR}/salmon/mRNA/${baseName} --gcBias --numBootstraps 100

# quantify ncRNA transcripts
${SALMON_SING} salmon quant -i ${ncRNAindex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz -o ${WORK_DIR}/salmon/ncRNA/${baseName} --gcBias --numBootstraps 100

# quantify pseudoRNA transcripts
${SALMON_SING} salmon quant -i ${pseudoIndex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz -o ${WORK_DIR}/salmon/pseudoRNA/${baseName} --gcBias --numBootstraps 100

# quantify TnRNA transcripts
${SALMON_SING} salmon quant -i ${tnIndex} -l A -r ${WORK_DIR}/cutadapt/${baseName}.fastq.gz -o ${WORK_DIR}/salmon/tnRNA/${baseName} --gcBias --numBootstraps 100


