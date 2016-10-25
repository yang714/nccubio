#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
install.packages("seqinr")
library("seqinr")
library(Biostrings)
#pam1<-read.table("PAM250.txt")
#pam2<-read.table("PAM100.txt")
#dim(pam1)
#str(pam1)
#dim(pam2)
#str(pam2)
##########
#s1 <- "ARR"
#s2 <- "RRR"
#########
#class(pam1)
seqs<-read.fasta("test.fasta" )
# construct PAM250 from PAM1
#pam250<-as.matrix(pam1)
#pam100<-as.matrix(pam2)

data(PAM250)
PAM250
dim(PAM250)
#data.matrix(PAM100)
#PAM100
#dim(PAM100)



#PAM250[1:2,1:3]
#test<-seqs[[2]]
#choosebank("swissprot")
lep1 <- getSequence(seqs[[1]])
lep1
lep2 <- getSequence(seqs[[2]])
lep2
lep1ing<-c2s(lep1)
lep2ing<-c2s(lep2)
lep1ing
lep1ing<- toupper(lep1ing)
lep1ing
lep2ing<- toupper(lep2ing)
lep2ing
#lepraeseq <- getSequence(seqs[[1]])
#lepraeseq
#s3 <- "GAT"
#s4 <- "GAT"
#sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
#almglobal <- pairwiseAlignment(lep1ing, lep2ing, type="global",substitutionMatrix = "PAM250", scoreOnly = FALSE)
almlocal <- pairwiseAlignment(lep1ing, lep2ing, type="local",substitutionMatrix = "PAM250", scoreOnly = FALSE)
almglobal <- pairwiseAlignment(lep1ing, lep2ing, type="global",substitutionMatrix="PAM250", scoreOnly = FALSE)
#alm <- pairwiseAlignment(lep1, lep2, scoreOnly = FALSE)
#alm <- pairwiseAlignment(s3, s4,type="global",substitutionMatrix = "PAM250", scoreOnly = FALSE)
#alm
#sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
#almlocal100 <- pairwiseAlignment(lep1ing, lep2ing, type="local",substitutionMatrix = "PAM100", scoreOnly = FALSE)
#almglobal100 <- pairwiseAlignment(lep1ing, lep2ing, type="global",substitutionMatrix="PAM100", scoreOnly = FALSE)
#alm <- pairwiseAlignment(s3, s4, type="global", substitutionMatrix =sigma,  scoreOnly = FALSE)
#alm <- pairwiseAlignment(s3, s4, type="global",  scoreOnly = FALSE)
almglobal
almlocal 
p1<-as.character(pattern(almglobal))
p2<-as.character(subject(almglobal))
p3<-c(p1,p2)

p4<-as.character(pattern(almlocal))
p5<-as.character(subject(almlocal))
p6<-c(p4,p5)

p7<-c(p3,p6)


#g1<-as.character(pattern(almglobal100))
#g2<-as.character(subject(almglobal100))
#g3<-c(g1,g2)

#g4<-as.character(pattern(almlocal100))
#g5<-as.character(subject(almlocal100))
#g6<-c(g4,g5)

#g7<-c(g3,g6)



#write.table(p3,"ppp.txt")
write.fasta(as.list(p7),names=names(p7),file.out="result250.fasta")
#write.fasta(as.list(g7),names=names(g7),file.out="result100.fasta")
#r<-as.character(almglobal)
#dna0 <- DNAStringSet
#name<-as.list(almglobal)
#write.fasta(sequences=almglobal,names=names(almglobal),file.out="write.my.dna.fasta")
#r=c(subject(almglobal),pattern(almglobal))
#r = c(as(subject(almglobal), "DNAStringSet"),  as(pattern(almglobal), "DNAStringSet"))
#names(r) = c("subject", "pattern")
#writeXStringSet(r,"out.fasta")ile
#write.(name,"rrrrrr.txt")
#write.fasta(as.list(name),names=names(name),file.out="rrrrrr.fasta")

