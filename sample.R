#Computation for HBV caccination data
#methylation
sample <- read.csv("sample.csv",sep="\t",header=F)
index <- order(as.numeric(gsub("_V","",gsub("GR","",sample[,2]))))
x_methyl <- read.csv("GSE161020_series_matrix.txt.gz",sep="\t",comment.char="!")
X <-  t(data.matrix(scale(x_methyl[,-1][,index]))) %*% data.matrix(scale(x_methyl[,-1][,index]))
#Gene expression
files <- list.files("./",pattern="count.txt.gz")
x_all <- NULL
for (i in c(1:length(files)))
{
    cat(i," ")
    x <- read.csv(files[i],sep="\t",header=F)
    x_all <- cbind(x_all,x[,2])
}
x_all <- data.frame(x[,1],x_all)
X0 <- t(data.matrix(scale(x_all[,-1]))) %*% data.matrix(scale(x_all[,-1]))
#proteome
x <- read.csv("GR01,04,09,10,11,13,15,17,18,19.txt",sep="\t")
x1 <- read.csv("GR02,03,05,06,07.txt",sep="\t")
x[1:3,2] <- c("a","b","c")
x1[1:3,2] <- c("a","b","c")
x0 <- merge(x,x1,by.x=2,by.y=2,all=T,sort=F)
cells <- unlist(lapply(strsplit(colnames(x0),".",fixed=T),"[",1))
GR <- unlist(x0[1,])
visit <- unlist(x0[2,])
labels_cells <- names(table(cells))[-3]
labels_GR<- names(table(GR))[-c(1:2)]
labels_visit <- names(table(visit))[-c(1:2)]
visit[!(visit %in% labels_visit)] <- "X"
cells[!(cells %in% labels_cells)] <- "X"
GR[!(GR %in% labels_GR)] <- "X"


x_all_WBC <- NULL
x_all_Plasma <- NULL

for (i in c(1:length(labels_GR)))
{
    cat(i," ")
    for (j in c(1:length(labels_visit)))
    {
        index_WBC <- visit==labels_visit[j] & GR==labels_GR[i] & cells=="WBC"
        x_all_WBC = cbind(x_all_WBC,x0[,index_WBC,drop=F][,1])
        index_plasma <- visit==labels_visit[j] & GR==labels_GR[i] & cells=="Plasma"
        x_all_Plasma = cbind(x_all_Plasma,x0[,index_plasma,drop=F][,1])
    }
}
x_all_WBC <- data.frame(x0[,1],x_all_WBC)
x_all_Plasma <- data.frame(x0[,1],x_all_Plasma)
x_all_WBC[is.na(x_all_WBC)]<-0
x_all_Plasma[is.na(x_all_Plasma)]<-0
x_all_WBC[,-1] <- apply(x_all_WBC[,-1],2,as.numeric)
x_all_WBC[is.na(x_all_WBC)]<-0
x_all_WBC[,-1] <- scale(x_all_WBC[,-1])
x_all_WBC[is.na(x_all_WBC)]<-0
X1_WBC <- t(data.matrix(x_all_WBC[,-1])) %*% data.matrix(x_all_WBC[,-1])
x_all_Plasma[,-1] <- apply(x_all_Plasma[,-1],2,as.numeric)
x_all_Plasma[is.na(x_all_Plasma)]<-0
x_all_Plasma[,-1] <- scale(x_all_Plasma[,-1])
x_all_Plasma[is.na(x_all_Plasma)]<-0
X1_Plasma <- t(data.matrix(x_all_Plasma[,-1])) %*% data.matrix(x_all_Plasma[,-1])
dim(X) <- c(5,15,75)
dim(X0) <- c(5,15,75)
dim(X1_WBC) <- c(5,15,75)
dim(X1_Plasma) <- c(5,15,75)
Z <- array(NA,c(75,75,4))
Z[,,1] <- X/mean(X)
Z[,,2] <- X0/mean(X0)
Z[,,3] <- X1_WBC/mean(X1_WBC)
Z[,,4] <- X1_Plasma/mean(X1_Plasma)
dim(Z) <- c(5,15,5,15,4)
require(rTensor)
HOSVD <- hosvd(as.tensor(Z))
u <- outer(HOSVD$U[[1]][,2],HOSVD$U[[2]][,1])
dim(u) <- 75
P <- pchisq(scale(data.matrix(scale(x_methyl[,-1][,index])) %*% u)^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
y <- read.csv("GPL21145_MethylationEPIC_15073387_v-1-0.csv.gz",sep=",",skip=7)
gene <- y[match(x[p.adjust(P,"BH")<0.01,1],y[,1]),15]
gene <- sort(unique(unlist(strsplit(gene,";"))))
P <- pchisq(scale(data.matrix(scale(x_all[,-1])) %*% u)^2,1,lower.tail=F)
gene <- x_all[p.adjust(P,"BH")<0.01,1] #genes selected by mnethylation
P <- pchisq(scale(data.matrix(x_all_WBC[,-1]) %*% u)^2,1,lower.tail=F)
gene <- x_all[p.adjust(P,"BH")<0.01,1] #genes selected by gene expression
gene <- unlist(strsplit(x_all_WBC[p.adjust(P,"BH")<0.05,1],"|",fixed=T)) 
gene <- gene[grep("NP",gene)] 
gene <- unlist(lapply(strsplit(gene,".",fixed=T),"[",1)) #genes selected by WBC
P <- pchisq(scale(data.matrix(x_all_Plasma[,-1]) %*% scale(u))^2,1,lower.tail=F)
gene <- unlist(strsplit(x_all_Plasma[p.adjust(P,"BH")<0.05,1],"|",fixed=T))
gene <- gene[grep("NP",gene)]
gene <- unlist(lapply(strsplit(gene,".",fixed=T),"[",1)) #genes selected by Plasma

#Synthetic data set
label <- matrix(F,ncol=3,nrow=1000)
  label[1:10,]<-T
  for (j in c(1:3))
  {
      label[10*j+c(1:10),j] <-T
  }
 P1 <- NULL
 require(rTensor)
 xx <- seq(0,1,length=10)
for (k in c(1:100))
{
    cat(k," ")
  x <- matrix(runif(1000*10),10,1000)
  x1 <- matrix(runif(1000*10),10,1000)
  x2 <- matrix(runif(1000*10),10,1000)
  x[,1:10] <- x[,1:10]+xx
  x1[,1:10] <- x1[,1:10]+xx
  x2[,1:10] <- x2[,1:10]+xx
  x[,11:20] <- x[,11:20]+xx
  x1[,21:30] <- x1[,21:30]+xx
  x2[,31:40] <- x2[,31:40]+xx
  x <- t(x)
  x1 <- t(x1)
  x2 <- t(x2)
  Z <- array(NA,c(1000,10,3))
  Z[,,1] <- x
  Z[,,2] <- x1
  Z[,,3] <- x2
  dim(Z) <- c(1000,30)
  ZZ <- t(Z) %*% Z
  dim(ZZ) <- c(10,3,10,3)
  HOSVD <- hosvd(as.tensor(ZZ))
  plot(HOSVD$U[[1]][,2])
  plot(HOSVD$U[[2]][,1])
  U <- Z %*% as.vector(outer(HOSVD$U[[1]][,2],HOSVD$U[[2]][,1]))
  plot(U)
  P <- pchisq(scale(U)^2,1,lower.tail=F)
  P1 <- cbind(P1,P)
  }
  
 PP1 <- apply(P1,2,function(x){p.adjust(x,"BH")})
 sum(rowSums(PP1<0.01)[1:10])
 
 sum(rowSums(PP1<0.01)[11:1000])
