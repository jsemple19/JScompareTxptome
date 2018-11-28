#! /usr/bin/make -f
# mapping Ringo's RNAseq data
#

###############################
########### VARIABLES #########
###############################

genomeFile=${HOME}/genomeVer/ws260/sequence/c_elegans.PRJNA13758.WS260.genomic.fa
genomeDir=$(dir ${genomeFile})
# annotFile:=  /home/ubelix/izb/semple/genomeVer/ws260/annotation/c_elegans.PRJNA13758.WS260.annotations.gff3.gz
# gunzip $annotFile
# annotFile=${annotFile%.gz}
# use gffread from cufflinks to convert gff to gtf
# module add UHTS/Assembler/cufflinks/2.2.1
# b=(`basename -s .gff3 ${annotFile}`)
# gffread $annotFile -T -o ${annotFile%/*}/${b}.gtf
# need to remove wierd exons:
# grep WormBase c_elegans.PRJNA13758.WS260.annotations.gtf > c_elegans.PRJNA13758.WS260.annotations1.gtf
# mv c_elegans.PRJNA13758.WS260.annotations1.gtf c_elegans.PRJNA13758.WS260.annotations		.gtf
annotFile=${HOME}/genomeVer/ws260/annotation/c_elegans.PRJNA13758.WS260.annotations.gtf
#mRNAseqFile=${genomeDir}/c_elegans.PRJNA13758.WS260.mRNA_transcripts.fa.gz
mRNAindex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_mRNA_index
ncRNAindex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_ncRNA_index
pseudoIndex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_pseudogenic_index
tnIndex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_transposon_index


trimmomaticDIR=/software/UHTS/Analysis/trimmomatic/0.36/bin
#trimAdapterFile=${HOME}/adaptors/TruSeq_3_SE.fa
trimAdapterFile=/software/UHTS/Analysis/trimmomatic/0.36/bin/adaptors/TruSeq3-SE.fa
starPath=/software/UHTS/Aligner/STAR/2.6.0c/bin/
#trimmomaticDIR := ${HOME}/Trimmomatic-0.36
#trimAdapterFile := ${trimmomaticDIR}/adapters/TruSeq_2-3_PE.fa
#bname :=  $(addprefix 180126_SNK268_A_L001_JIB-, 1 2 3 4)
#longbname := $(addsuffix _R1, $(bname)) $(addsuffix _R2, $(bname))
#PRESEQ:=/Users/semple/mySoftware/preseq
#rawDirName := /home/ubelix/izb/semple/Francesca/Ringo_RNA_seq/rawData
#fileListFile=../fastqList.txt
#fileList= ( $(cat ../fastqList.txt) )


parentDir=$(dir ${fastqFile})
analysisDir=$(parentDir:%rawData/=%JScompareTxptome)

fileName:=$(notdir ${fastqFile})
baseName:=$(basename $(basename ${fileName})) #to remove both .fastq and .gz
#baseName= `basename -s .fastq.gz -a ${fastqFile}`

#newList=${newList[@]/.fastq.gz/.bam}


#list of the final output files
objects :=$(addsuffix _fastqc.html, ../rawData/fastQC/${baseName})\
	$(addsuffix _fastqc.html, cutadapt/fastQC/${baseName})\
	$(addsuffix _fastqc.html, trim/fastQC/${baseName})\
	$(addsuffix .Aligned.out.bam, bam/${baseName}) \
	$(addprefix salmon/mRNA/,$(addsuffix /quant.sf,${baseName})) \
	$(addprefix salmon/ncRNA/,$(addsuffix /quant.sf,${baseName})) \
	$(addprefix salmon/pseudoRNA/,$(addsuffix /quant.sf,${baseName})) \
	$(addprefix salmon/tnRNA/,$(addsuffix /quant.sf,${baseName})) 

intermediateObjects:=$(addprefix cutadapt/, $(addsuffix .fastq.gz, ${baseName}))\
	$(addprefix trim/, $(addsuffix .fastq.gz, ${baseName})) \
	$(addsuffix .Aligned.out.sam, bam/${baseName}) \

###############################
########### RULES  ############
###############################


all: $(objects) $(intermediateObjects)

print-%  : ; @echo $* = $($*)

#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on sequences
../rawData/fastQC/%_fastqc.html: ../rawData/%.fastq.gz
	mkdir -p ../rawData/fastQC
	fastqc $^ -o ../rawData/fastQC 
	
#######################################################
## trim adaptors with cutadapt                       ##
#######################################################

# use cutadapt to trim
cutadapt/%_R1.fastq.gz: ../rawData/%_R1.fastq.gz
	mkdir -p cutadapt
	#adaptor=`grep $*_R1.fastq.gz ../sampleList.csv | cut -d"," -f7
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
		-o cutadapt/$*_R1.fastq.gz  \
		../rawData/$*_R1.fastq.gz

#redo fastQC on trimmed reads
cutadapt/fastQC/%_fastqc.html: cutadapt/%.fastq.gz
	mkdir -p cutadapt/fastQC
	fastqc cutadapt/$*.fastq.gz -o cutadapt/fastQC 

#######################################################
## quality trim reads with Trimmomatic               ##
#######################################################

# there is evidence that quality trimming might make RNAseq alignments worse
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2
# Trimming of sequence reads alters RNA-Seq gene expression estimates
# CA Williams et al.(2016)
# But the trimmomatic trimming i use is not very aggressive 

# use trimmomatic to trim
trim/%_R1.fastq.gz: cutadapt/%_R1.fastq.gz
	mkdir -p trim
	java -jar ${trimmomaticDIR}/trimmomatic-0.36.jar SE cutadapt/$*_R1.fastq.gz trim/$*_R1.fastq.gz ILLUMINACLIP:${trimAdapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> trim/report_$*_trimmomatic.txt

# redo fastQC on trimmed reads	
trim/fastQC/%_fastqc.html: trim/%.fastq.gz
	mkdir -p trim/fastQC
	fastqc trim/$*.fastq.gz -o trim/fastQC
		
#######################################################
## Align to genome with STAR                         ##
#######################################################


# index genome
# STAR --runMode genomeGenerate --genomeDir ${genomeFile%/*} --genomeFastaFiles ${genomeFile} --sjdbGTFfile ${annotFile} --runThreadN 4

# align to genome
bam/%.sam: trim/%.fastq.gz
	mkdir -p bam
	STAR --genomeDir ${genomeDir}  --readFilesIn trim/$*.fastq.gz --readFilesCommand zcat --outFileNamePrefix bam/$*. --runThreadN 4 --alignIntronMax 500 --quantMode GeneCounts

# convet sam to bam
bam/%.bam: bam/%.sam
	samtools view -b $^ -o $@
	rm $^

#######################################################
## Count reads with Salmon                           ##
#######################################################

# prepare reference genome
# See separate script indexTxpts4salmon.sh

######## NOTE: I do not have estimates for --fldMean and --fldSD as i have no access to the bioanalyser files #########
#######  therefore the quantification based on the effective transcript length will be wrong!!!! ######

# quantify mRNA transcripts
salmon/mRNA/%/quant.sf: trim/%.fastq.gz
	salmon quant -i ${mRNAindex} -l A -r trim/$*.fastq.gz -o salmon/mRNA/$* --gcBias --numBootstraps 100

# quantify ncRNA transcripts
salmon/ncRNA/%/quant.sf: trim/%.fastq.gz
	salmon quant -i ${ncRNAindex} -l A -r trim/$*.fastq.gz -o salmon/ncRNA/$* --gcBias --numBootstraps 100

# quantify pseudoRNA transcripts
salmon/pseudoRNA/%/quant.sf: trim/%.fastq.gz
	salmon quant -i ${pseudoIndex} -l A -r trim/$*.fastq.gz -o salmon/pseudoRNA/$* --gcBias --numBootstraps 100

# quantify TnRNA transcripts
salmon/tnRNA/%/quant.sf: trim/%.fastq.gz
	salmon quant -i ${tnIndex} -l A -r trim/$*.fastq.gz -o salmon/tnRNA/$* --gcBias --numBootstraps 100


