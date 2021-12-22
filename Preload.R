# Augmat <- function(data.model, T.max , time.adjust= T) 
#          {
#          X0<-diag(T.max)
#          n <- nrow(data.model)
#          X <- data.model[,-c(1:2)]
#          p <- ncol(data.model)-2
#          X.temp<-list()
#          for(i in 1:n)
#          X.temp[[i]]<-cbind(X0,t(matrix(ncol=T.max,nrow=p,rep(X[i,],T.max))))
#          
#          X.final<-numeric()
#          y<-numeric()
#          id<-numeric()
#          time<-numeric()
#           
#          data.model$event <-data.model[,1]
#          if(time.adjust){data.model$event <- data.model$event+1}
#            
#          for(i in 1:n){
#            y.i <- c(rep(0,(data.model$event[i]-1)),data.model[i,2])
#            time.i<-1:data.model$event[i]
#            X.final<-rbind(X.final,X.temp[[i]][time.i,])
#            y<-c(y,y.i)
#            id<-c(id,rep(i,data.model$event[i]))
#            time<-c(time,time.i)
#          }
#          # The following step is necessary to format the every column as numeric
#           X.mat = matrix(NA, nrow=nrow(X.final), ncol=ncol(X.final))
#           for(i in 1:ncol(X.final)){X.mat[,i] <- as.numeric(X.final[,i])}
#          data.Aug<-cbind(data.frame(y=y,time=time, id=id),X=X.mat)
#          data.Aug$id<-as.factor(data.Aug$id)
#          table(table(data.Aug$id))
#          return(data.Aug)
# }


GenData <- function(minfreq=3, maxfreq = 29,n=100,p=500, Left.trun.prob=NULL,rho.vec= NULL,NZ= c(-4,-4,-4,8,8),nCat=0,... ){ 
test.ind<-T
k<- 1;
while(test.ind & k<=20)
{
  Simout<- Sim_fun( n=n,p=p,rho.vec= rho.vec, NZ= NZ,Left.trun.prob=Left.trun.prob, ...)
  
  data_s<- Simout$data_s
  
  if( all(table(data_s$time) >  (minfreq-1))  & all(table(data_s$time)< (maxfreq+1)))
  {
    test.ind<-F
    print (paste0("took ", k , " many steps"))
  }
  restr<- as.vector(table(data_s$time))
  if(min(restr) < minfreq)
  {
    print(paste0("Min. frequency  <  ", minfreq))
    # print(restr) 
  }else{if(max(restr) > maxfreq)print(paste0("Max. frequency > ", maxfreq))}
  k <- k+1
}
Simout
}

Sim_fun= function(n,p,rho.vec= c(0.7,0.4),baseline=NULL,FxdRgtCen=F,block=c(3,2), link="logit",  Left.trun.prob=c(0.6,0.2,0.2), silent=F, cencor.prob= 0.05, NZ= c(-4,-4,-4,8,8), seed=NULL,frailty_Q=1, nonL=F ,beta_cat=NULL, T.max=10)
{
  
  library(MASS)
  library(Matrix)
  ################################
  
  
  family=binomial(link=link)
  
  if(!is.null(seed))set.seed(seed)
  
  b<-rnorm(n,0,sd=frailty_Q)  #random effect
  
  
  
  test.seq<-seq(1,T.max,by=1)
  #baseline<-rep(-2,10)
  if(is.null(baseline)){baseline<-dgamma(test.seq-2, shape=5, scale = 1)*2-2.3}
  
  # Using the normal vector
  # X<-matrix(ncol=p,nrow=n,runif(p*n,min=0,max=1))
  #   sigma1 <- t(chol(0.9*matrix(1, nrow=3,ncol=3)+(1-0.9)*diag(3)))
  #   sigma2 <- t(chol(0.75*matrix(1, nrow=2,ncol=2)+(1-0.75)*diag(2)))
  #   sigmaL  <- bdiag(sigma1, sigma2)
  #   normat <- (matrix(ncol=5,nrow=n, rnorm( n *5, mean=-1, sd=1 )) %*% t(sigmaL))
  #  X[,1:5] <- as.matrix(normat)
  #  
  # Using the uniform vector
  if(is.null(rho.vec)){
    X<-matrix(ncol=p,nrow=n,runif(p*n,min=0,max=1))
  }else{
    nblk <- length(rho.vec)
    X1<-matrix(ncol=(p+nblk),nrow=n,runif((p+nblk)*n,min=0,max=1))
    if(nblk==1){block= sum(block)}
    X=X1[,(nblk+1):(p+nblk)]
    k1<-1
    for(jj in 1:nblk){
      size <- block[jj]
      X[,k1:(k1+size-1)]<-block.corr(rho.vec[jj],X1[,c(jj,(k1+nblk):(k1+nblk+size-1))], size)
      k1<- k1+size
    }
  }
  # cat("testing1\n")
  #eta.test<-X%*%c(3,5,-10,rep(0,p-3))
  #summary(eta.test)
  
  X0<-diag(T.max)
  
  if(nonL){
    if(is.null(beta_cat)){beta_cat <- c(5,-5,6,-4);}
    Des <- polyDes(X[,p],3)
    more<-list(spl= X[,p]);
    more$cor <- cor(Des[,-1]);
    Des1 <- bs.design(X[,p], diff.ord=2, spline.degree=3, knots.no=10)$B
    beta<-c(baseline,beta_cat,NZ,rep(0,p-(length(NZ)+length(beta_cat))))
    X <- cbind(Des,X[,1:(p-length(beta_cat))])
    #print(dim(X))
  }else{
    beta<-c(baseline,NZ,rep(0,p-length(NZ)))
    more=NA;
  }
  
  X.temp<-list()
  for(i in 1:n)
    X.temp[[i]]<-cbind(X0,t(matrix(ncol=T.max,nrow=p,rep(X[i,],T.max))))  
  
  X.final<-numeric()
  
  y<-numeric()
  id<-numeric()
  time<-numeric()
  L<-numeric()
  eta.vec<-c()
  Res <- c()
  
  
  i<-1
  while(i  <= n)
  {
    mu.i<-numeric() 
    k<-1
    next.try<-T
    y.i<-numeric()
    
    while(k<T.max+1 && next.try){
      next.try<-F
      eta<-X.temp[[i]][k,]%*%beta+b[i]
      eta.vec[i]<-eta
      mu.i<-c(mu.i,family$linkinv(eta))
      p.hit<-mu.i[k]#prob(mu.i)
      y.temp<-rbinom(1,1,p.hit)
      
      if(FxdRgtCen){   
        if(y.temp!=1)next.try<-as.logical(1)
      }else{
        if(y.temp!=1)next.try<-as.logical(1-rbinom(1,1,cencor.prob))
      }
      y.i<-c(y.i,y.temp)
      k<-k+1
    }
    li=1;
    if(!is.null(Left.trun.prob)){
      li <-sample(x=length(Left.trun.prob), size=1, replace=TRUE, prob=Left.trun.prob)
      #print(c(w.i,(k-1)))
      if(li < k){
        time.i <- li:(k-1)
        X.final<-rbind(X.final,X.temp[[i]][li:(k-1),])
        time<-c(time,time.i)
        delta.i <- y.i[k-1] 
        y.i<-y.i[li:(k-1)]
        y<-c(y,y.i)
        id<-c(id,rep(i,(k-li)))
        Res <- rbind(Res, c(i,k-1,delta.i,li))
        i<-(i+1)
      }else{ if(!silent)cat(paste0(i,"-th simulation: rejected\n "))}
    }else
    {               
      time.i<-1:(k-1)
      X.final<-rbind(X.final,X.temp[[i]][li:(k-1),])
      y<-c(y,y.i)
      delta.i <- y.i[k-1]
      id<-c(id,rep(i,k-1))
      time<-c(time,time.i)
      Res <- rbind(Res, c(i, k-1,delta.i,li))
      i<-i+1
    }
    
    
  }
  data_long<-data.frame(id, y,time=time,X=X.final)
  data_short<-data.frame(Res,X=X)
  colnames(data_short)[1:4]  <- c("id","time", "delta", "Lt")
  
  data_long$id<-as.factor(data_long$id)
  data_short$id<-as.factor(data_short$id)
  
  out=list(family=family,data_s=data_short,data_l= data_long, baseline=baseline, beta=beta, test.seq=test.seq, eta=eta.vec, more=more)
}


WideData = function(data)
{ # Uses the output data from Simdata:
  T.max  <- max(data$time)
  X0<-diag(T.max)
  X  <- as.matrix(data[, -(1:4)])
  p <- ncol(X)
  X.final<-numeric()
  Lt <- data$Lt
  id <- c()
  y <-c()
  time <- c() 
  for(i in 1:length(data$id)){
    id <- c(id,rep(data$id[i],(data$time[i]-data$Lt[i]+1)))
    
    if((data$time[i]-data$Lt[i]) >0){y.i <- c(rep(0,data$time[i]-data$Lt[i]), data$delta[i])}else{y.i <-  data$delta[i]}
    y <- c(y,y.i);
    time <- c(time, data$Lt[i]:data$time[i]) 
    X.temp<-cbind(X0,t(matrix(ncol=T.max,nrow=p,rep(X[i,],T.max))))
    X.final<-rbind(X.final,X.temp[data$Lt[i]:data$time[i],])
  }
  data_long<-data.frame(id, y,time=time,X=X.final)
  data_long$id<-as.factor(data_long$id)
  data_long
}


sourceDir <- function(path, trace = F, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

CrossParamTab <- function(parnames,lst){ 
  require(tidyverse); 
  for(i in 1:length(parnames)){
    if(i >1){l1= nrow(tt) ; l2= length(lst[[parnames[i]]]); 
    tt <-cbind.data.frame( tt[rep(1:l1,each=l2),], rep(lst[[parnames[i]]],l1), stringsAsFactors = F)
    colnames(tt)[i] <- parnames[i]}else{
      tt <- data.frame(X1 =lst[[parnames[1]]], stringsAsFactors = F)
    }}
  colnames(tt)[1] <- parnames[1]
  tt}

block.corr<- function(rho, mat, k=2)
{
  if( ncol(mat) < k+1 )stop(paste0("number of columns must be strictly greater than ", k))
  rho.1<- sqrt(rho)
  beta <- rho.1/sqrt(1-rho)
  out= mat[,2:(k+1)]+ beta*mat[,1]
  return(out)
}


 