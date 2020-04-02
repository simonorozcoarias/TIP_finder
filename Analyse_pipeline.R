#####
#Date : 20 décembre 2016
# Analyse automatique des résultats du pipeline 3000genomes
#####

#elimination des warning
options(warn=-1)

#recuperation des arguments
args <- commandArgs(trailingOnly = TRUE)

#TE family
te<-as.character(args[1])

inputFile<-try(system("ls *txt",intern=TRUE))
files <- unlist(strsplit(inputFile,split=" "))

length(files)
M<-data.frame()

all<-paste("all_insertion_",te,".names",sep="")
M<-read.table(all,sep="\n",h=F)

for (i in 1:length(files)){
 res<-read.table(files[i],sep="\n",h=F)
 tmp<-M$V1 %in% res$V1
 M<-cbind(M,tmp)
 }
colnames(M)<-c("insertion",files)

#pos cov5
pos<-paste("all_position_cov5_",te,".names",sep="")
POS<-read.table(pos,sep="\n",h=F)

#Matrice final
M2<-M[which(M$insertion %in% POS$V1),]
write.table(M2,file="matrice_final.csv",sep="\t",row.names=F,col.names=T,quote=F)



