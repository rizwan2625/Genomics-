##PART I###########USE FUNCTION EIGEN##########################
library(matlib)

#Generate a matrix with 6 genes and 6 samples containing data from a normal distribution 
names<-list(c("gene1","gene2","gene3","gene4","gene5","gene6"),c("sample1","sample2","sample3","sample4","sample5","sample6"))
A<-matrix(rnorm(36),nrow=6,ncol=6,dimnames=names) #do not re-run A for part II, so you can compare part I and II
A

#Calculate the covariance matrix
#S<-scale(A,center=TRUE,scale=FALSE) #By default, this function will standardize the data (mean zero, unit variance). To indicate that we just want to subtract the mean, we need to turn off the argument scale = FALSE.
C<-cov(A)#sames as C<-cov(S)
det(C)# so it is diagonalizable
E<-eigen(C)

#Calculate the eigenvalues and proportion of variance
round(inv(E$vectors)%*%C%*%E$vectors,digits=2)
round(E$values,digits=2)
weight<-((E$values/sum(E$values))*100)
plot(weight, type="s", ylab="Percentage of variance", xlab="PCA component")

#Calculate the eigen vectors and plot transformed dataset
loadings<-E$vectors
names2<-list(c("sample1","sample2","sample3","sample4","sample5","sample6"),c("PC1","PC2","PC3","PC4","PC5","PC6"))
dimnames(loadings)<-names2
round(loadings,digits=3)
plot(loadings)
text(loadings,labels=row.names(loadings),pos=c(2,2,2,4,4,2))#please arrange labels so you can see sample names

#Calcualte the scores and plot them
S<-scale(A,center=TRUE,scale=FALSE)
Scores<-S%*%loadings
Scores
plot(Scores)
text(Scores,labels=row.names(Scores),pos=c(2,2,2,2,4,2))#please arrange labels so you can see gene names names

##PARTII###############USE FUNCTION PRINCOMP###########################

#calculate PCA
PCA<-princomp(A,cor=FALSE,scores=TRUE)#why do we set cor=FALSE?
PCA[1:7]#make sure you know what all data means, 
#tip, for PCA$center compare it to S from part I

plot(PCA)
#Check you would know how to calculate percentage of variance
#variance<-((PCA$sdev^2/sum(PCA$sdev^2))*100)
#plot(variance, type="s", ylab="Percentage of variance", xlab="PCA component")
#does it look like the one in part I?

PCA$loadings #are the same values as in P?
#par(mfrow = c(1, 2))#we define a plot with 1 row and 2 slots
plot(PCA$loadings)
text(PCA$loadings,labels=row.names(PCA$loadings),pos=c(4,4,4,2,2,4))#please arrange labels so you can see sample names
#run "par"  and plot together with eigenvectors loadings calculated in part I
#Do they look alike? What are the possible differences due to?

PCA$scores
plot(PCA$scores)
text(PCA$scores,labels=row.names(PCA$scores),pos=c(4,4,4,2,2,4))#please arrange labels so you can see sample names
#plot together with scores calculated in part I
#plot loading and scores together. What are the most important genes contributing to sample separation?


