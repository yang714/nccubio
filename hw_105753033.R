# read PAM1 from data
pam1<-read.table("pam1.txt")

# check PAM1 data
dim(pam1)
str(pam1)
#class(pam1)

# construct PAM250 from PAM1
pam2<-as.matrix(pam1)
pam2=pam2/10000
#pam2<-round(pam2, digit = 3)
pam250<-pam2
#class(pam2)
for(ggc in c(1:249))
{
  pam250 <- pam250%*%pam2
  #Reduce("%*%", pam2, accumulate = TRUE)
  #pam2<-pam2+pam2
}
pam250<-round(pam250, digit = 2)
data<-pam250*100

#
# output PAM250 as a file
write.table(data,file="pam250.txt")

