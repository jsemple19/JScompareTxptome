# script to obtiain txptome fasta files from wormbase and index them for salmon
genomeDir=/home/ubelix/izb/semple/genomeVer/ws260/sequence/

curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS260.mRNA_transcripts.fa.gz -o ${genomeDir}/c_elegans.PRJNA13758.WS260.mRNA_transcripts.fa.gz

curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS260.ncRNA_transcripts.fa.gz -o ${genomeDir}/c_elegans.PRJNA13758.WS260.ncRNA_transcripts.fa.gz

curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS260.pseudogenic_transcripts.fa.gz -o ${genomeDir}/c_elegans.PRJNA13758.WS260.pseudogenic_transcripts.fa.gz

curl ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/PRJNA13758/sequence/transcripts/c_elegans.PRJNA13758.WS260.transposon_transcripts.fa.gz -o ${genomeDir}/c_elegans.PRJNA13758.WS260.transposon_transcripts.fa.gz

module add vital-it
module add UHTS/Analysis/salmon/0.11.2;
cd $genomeDir

salmon index -t ${genomeDir}/c_elegans.PRJNA13758.WS260.mRNA_transcripts.fa.gz -i ws260_mRNA_index

salmon index -t ${genomeDir}/c_elegans.PRJNA13758.WS260.ncRNA_transcripts.fa.gz -i ws260_ncRNA_index

salmon index -t ${genomeDir}/c_elegans.PRJNA13758.WS260.pseudogenic_transcripts.fa.gz -i ws260_pseudogenic_index

salmon index -t ${genomeDir}/c_elegans.PRJNA13758.WS260.transposon_transcripts.fa.gz -i ws260_transposon_index

# use these addresses to reference indeces in main script
mRNAindex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_mRNA_index
ncRNAindex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_ncRNA_index
pseudoIndex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_pseudogenic_index
tnIndex=/home/ubelix/izb/semple/genomeVer/ws260/sequence/ws260_transposon_index

