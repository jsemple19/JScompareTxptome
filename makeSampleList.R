
#l<-list.files("../rawData",pattern="fastq.gz",all.files=T)

sampleList<-read.delim("./ringoSampleList.txt"," ",stringsAsFactors=F,header=T)

newList<-as.data.frame(sampleList$Strain)
names(newList)<-"strain"
idx<-sampleList$HS=="HS"
newList$hs<-as.numeric(idx)
newList$fed<-FALSE+(sampleList$Fed=="fed")
newList$sampleID<-sampleList[,"Library_ID"]
newList$libType<-"stranded_mRNA"
newList$index<-sampleList[,"Index"]

write.csv(newList,"../libList.csv",row.names=F)

strainList<-data.frame(strain=c("200", "218", "284", "285"),
                       hsHLH1=c(1,0,0,1),
                       mes2=c(0,0,1,1))

write.csv(strainList,"../strainList.csv",row.names=F)


fastqList<-read.delim("../fastqList.txt",stringsAsFactors=F,header=F)
names(fastqList)<-"fastqFile"

nameFields<-sapply(basename(fastqList$fastqFile),strsplit,split="_")
fastqList$date<-sapply(nameFields,"[[",1)
fastqList$lane<-sapply(nameFields,"[[",4)
fastqList$sampleID<-sapply(nameFields,"[[",5)

idx<-match(fastqList$sampleID,newList$sampleID)
fastqList[,c("strain","hs","fed","index")]<-newList[idx,c("strain","hs","fed","index")]

idx<-match(fastqList$strain,strainList$strain)
fastqList[,c("hsHLH1","mes2")]<-strainList[idx,c("hsHLH1","mes2")]

write.csv(fastqList,"../sampleList.csv",row.names=F)
