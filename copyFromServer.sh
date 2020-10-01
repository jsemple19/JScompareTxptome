#! /bin/bash
# copy results from server

sourceDir=$1
destDir=$2
serverAddress=semple@submit.unibe.ch:

dirs2copy=(bam salmon ref makefile *.sh)

for d in ${dirs2copy[@]}
do
	scp -r ${serverAddress}${sourceDir}$d $destDir/
done
