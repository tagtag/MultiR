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
