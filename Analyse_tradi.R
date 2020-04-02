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

mat<-paste("matrice_final_",te,".csv",sep="")
M<-read.table(mat,sep="\t",h=T)
tmat<-as.data.frame(t(M))
tmat$varieties<-row.names(tmat)
REF<-tmat

tradi<-read.table("/home/mchristine/Bureau/3000g/all_names_3K_Traditional.txt",sep="\n",h=F)

REFT<-REF[which(REF$varieties %in% tradi$V1),]
colnames(REFT)<-c(as.vector(M[,1]),"varieties")
matT<-paste("matrice_final_tradi_",te,".csv",sep="")
write.table(REFT,file=matT,row.names=F,col.names=T,quote=F,sep="\t")

REFT<-read.table(matT,sep="\t",h=T)

nbv<-as.numeric(dim(REFT)[2])-1

RES<-colSums(REFT[,1:nbv])
png(paste("Nb_varieties_tradi_per_insertion_",te,".png",sep=""))
hist(log2(RES),nclass=50)
dev.off()

table(RES>1)

RES2<-rowSums(REFT[,1:nbv])
png(paste("Nb_insertion_per_varieties_tradi_",te,".png",sep=""))
hist(RES2,nclass=50)
dev.off()
