install.packages(seqinr)

library(seqinr)
library(Biostrings)
#x <- s2c("GAATTC")
#y <- s2c("GATTA")
#x <-s1<-"GAAT

#y <- s2<-"GATTA"


seqs<-read.fasta("test.fasta" )
#y <- s2c("GATTAG")
#x <- s2c("GAATTC") 
data(PAM250)

x <- getSequence(seqs[[1]])
x<-toupper(x)
x
y <- getSequence(seqs[[2]])
y<-toupper(y)
y
xx<-c2s(x)
yy<-c2s(y)
xx
yy
# calculate scoring matrix for matches (note the power of indexing)
#s <- BLOSUM62[y,x]
s<-PAM250[y,x]
s
#dim(s)
# set the gap penalty
d <-8
#cc<-PAM250
#cc
# compute max alignment score using dynamic programming

F <- matrix(data=NA,nrow=(length(y)+1),ncol=(length(x)+1))
############
PAM250[,2]
dim(PAM250)
PAM250$A

############

m<-length(y)*d
n<-length(x)*d
F[1,] <- -seq(0,n,d); F[,1] <- -seq(0,m,d)
rownames(F) <- c("",y); colnames(F) <- c("",x)
for (i in 2:(nrow(F)))
  for (j in 2:(ncol(F)))
  {
    F[i,j] <- max(c(F[i-1,j-1]+s[i-1,j-1],F[i-1,j]-d,F[i,j-1]-d))
    ################
   
    ################
  }
F[1,2]
FF<-F
F
dd1<-0
dd2<-0
dd3<-0
ccg<-0
#################

#################
almlocal <- pairwiseAlignment(xx, yy, type="local",substitutionMatrix = "PAM250", scoreOnly = FALSE)
almglobal <- pairwiseAlignment(xx, yy, type="global",substitutionMatrix="PAM250", scoreOnly = FALSE)
#almlocal
almglobal
###########################

p1<-as.character(pattern(almglobal))
p2<-as.character(subject(almglobal))
p3<-c(p1,p2)

p4<-as.character(pattern(almlocal))
p5<-as.character(subject(almlocal))
p6<-c(p4,p5)
FF<-as.matrix(F)
p7<-c(p3,p6)
#write.fasta(as.list(p7),names=names(p7),file.out="result250.fasta")

##################TRACEBACK

##################################

###################################################rtraceback-test
xo<-""
xo<- getSequence(xo[[1]])
yo<-""
yo<- getSequence(yo[[1]])
bad<-0
####################################################
yyy<-y
xxx<-x
see<-TRUE
sorcesee<-0
yyy[2]#G
F[2,2]
j<-2
i<-2
kx<-2
ky<-2
ppc<-0
te<-max(-4,-5,-9)
sorce<-0
maxl<-max(length(y),length(x))
for(k in 2:maxl)
{
  ppc<-ppc+1

  if (see==FALSE)
  {
    i<- i+1
    j<- j+1
    
  }
  
  
    else if(see==TRUE)
  {
   
  i<-i
  j<-j
    
  }
  

  
  
  if(max(c(F[i,j],F[i-1,j],F[i,j-1]))==F[i,j])
  {
  
    xo[k-1]<-xxx[kx-1]
    yo[k-1]<-yyy[ky-1]
    
    if(sorcesee==0)
      {
    sorce<-sorce+s[kx-1,ky-1]
    }
    else if(sorcesee==1){
      
      sorce<-sorce+s[kx-2,ky]
    }
    
    else if(sorcesee==2){
      
      sorce<-sorce+s[kx,ky-2]
    }
    
    
    i<-i+1
     j<-j+1
     kx<-kx+1
    ky<-ky+1
    
  
   
    see<-TRUE
    sorcesee==0
   
    
  }
  else if(max(c(F[i,j],F[i-1,j],F[i,j-1]))==F[i-1,j])
  {
  
     
    yo[k-1]<-"-"
    xo[k-1]<-xxx[kx-1]
    i<-i-1
    kx<-kx+1
    sorce<-sorce-8
    see<-FALSE
    sorcesee<-1
   
  }
  
  else if(max(c(F[i,j],F[i-1,j],F[i,j-1]))==F[i,j-1])
  {
    
  
####下
    
    xo[k-1]<-"-"
    yo[k-1]<-yyy[ky-1]
    j<-j-1
    sorce<-sorce-8
    
    ky<-ky+1
   see<-FALSE
   
   sorcesee<-2
  }
  
  else
  {
    bad+1
  }
  

}

#############################
#for(p in length(x):2)
 # xxx[p+1]<-xxx[p]
##############################
i
j
ppc
xo
length(xo)
yo

length(yo)
sorce
bad
#yyy
#xxx
see
xxo<-c2s(xo)
yyo<-c2s(yo)

xxo
yyo
p7<-c(xxo,yyo,sorce)
write.fasta(as.list(p7),names=names(p7),file.out="result.fasta")
#########################################################

###################
#for (i in 2:(nrow(F)))



