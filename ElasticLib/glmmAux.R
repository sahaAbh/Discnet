blockstand <- function(x, ipen.which, inotpen.which)
{
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50
  
  n <- nrow(x)
  x.ort <- x
  scale.pen <- list(); length(scale.pen) <- length(ipen.which)
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0){
    one <- rep(1, n)
    scale.notpen <- sqrt(drop(one %*% (x[,inotpen.which]^2)) / n)
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
  
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
           " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

correct.cat <- function(aaa,block)
{
  for(i in 1:length(block))
  {
    if(block[i]>1)
    {
      if(sum(aaa[(sum(block[1:(i-1)])+1):sum(block[1:i])])>0)
        aaa[(sum(block[1:(i-1)])+1):sum(block[1:i])] <- TRUE
    }  
  }
  return(aaa)
}


taylor.opt<-function(t_opt,y,X,fixef,ranef,Grad,family,P, K=NULL)
{
  delta<-c(fixef,ranef)+t_opt*Grad
  Eta<- X%*%delta
  if(is.null(K)){
    mu<-family$linkinv(Eta)
  }else{
    Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
    mu <- c(t(family$linkinv(Eta_cat)))
  }
  
  loglik <- logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE,K=K) 
  
  loglik<- -loglik + 0.5*t(delta[(length(fixef)+1):length(delta)]) %*% (P %*% delta[(length(fixef)+1):length(delta)])
  return(loglik)
}

taylor.opt.noRE<-function(t_opt,y,X,fixef,Grad,family, K=NULL)
{
  
  delta<-fixef+t_opt*Grad
  Eta<- X%*%delta
  if(is.null(K)){
    mu<-family$linkinv(Eta)
  }else{
    Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
    mu <- c(t(family$linkinv(Eta_cat)))
  }
  loglik <- -logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE, K = K) 
  
  return(loglik)
}


t.change<-function(grad,b)
{
  a<-(sign(b)==-sign(grad)&sign(grad)!=0)
  rate <- rep(Inf, length(a))
  rate[a] <--b[a]/grad[a]
  ret.obj<-list()
  ret.obj$min.rate<-min(rate)
  ret.obj$whichmin <- which.min(rate)
  return(ret.obj)
}

## Additions ##
gradient.elastic<-function(score.beta,b,lambda.b, alpha.b= alpha)
{
  p<-length(b)
  grad.elastic<-rep(0,p)
  # testing b[1:2]=c(1,2)
  b.isnt.0<-b!=0
  grad.elastic[b.isnt.0]<-score.beta[b.isnt.0]-lambda.b*(1-alpha.b)*b[b.isnt.0]-lambda.b*alpha.b*sign(b[b.isnt.0])
  b.is.0<-(b==0 & abs(score.beta)>= (lambda.b*alpha.b))
  grad.elastic[b.is.0]<-score.beta[b.is.0]-lambda.b*alpha.b*sign(score.beta[b.is.0])
  b.is.0<-(b==0 & abs(score.beta)<(lambda.b*alpha.b))
  grad.elastic[b.is.0]<-0
  return(grad.elastic)
}

gradient.elastic.block<-function(score.beta,b,lambda.b,alpha.b,block)       
{
  p<-length(b)
  grad.elastic<-rep(0,p)
  
  block2 <- rep(block,block)
  lambda_vec<-rep(lambda.b,length(b))
  lambda_vec <- lambda_vec*sqrt(block2)
  
  group.sum<-rep(0,length(block))
  group.sum[1]<-sqrt(sum(b[1:block[1]]^2))
  if(length(block)>1)
  {  
    for (i in 2:length(block))
      group.sum[i]<-sqrt(sum(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))
  }
  if(group.sum[1]!=0)
  {
    grad.elastic[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*((alpha.b*b[1:block[1]]/group.sum[1])+((1-alpha.b)*b[1:block[1]]))
  }else{
    
    if(group.sum[1]==0 & sqrt(sum(score.beta[1:block[1]]^2))> (alpha.b*lambda_vec[1]) )
    {
      grad.elastic[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*(alpha.b*score.beta[1:block[1]]/sqrt(sum(score.beta[1:block[1]]^2)))
    }else{
      grad.elastic[1:block[1]]<-0
    }}
  
  if(length(block)>1)
  {  
    for (i in 2:length(block))
    {
      if(group.sum[i]!=0)
      {
        grad.elastic[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*((alpha.b*b[(sum(block[1:(i-1)])+1):sum(block[1:i])]/group.sum[i])+ ((1-alpha.b)*b[(sum(block[1:(i-1)])+1):sum(block[1:i])]))
      }else{
        
        if(group.sum[i]==0 & sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))> (alpha.b*lambda_vec[sum(block[1:(i-1)])+1]))
        {
          grad.elastic[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*(alpha.b*score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]/sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2)))
        }else{
          grad.elastic[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-0
        }}
    }}
  return(grad.elastic)
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


summary.Discnet <- function(object, ...)
{
  se <- object$StdError
  zval <- coefficients(object) / se
  TAB <- cbind(Estimate = coefficients(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  res <- list(coefficients=TAB,baseline.eff= object$baseline ,Q=object$Q[2,])
  class(res) <- "summary.Discnet"
  res
}

