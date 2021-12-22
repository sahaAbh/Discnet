Discnet <- function(data_y,data_x, al = c(1), nus =25, index=NULL, Measure=c("bic", "perm"),lambda.large=NULL, init.theta=NULL, init.baseline=NULL,lambda.perm.length=500, family=binomial(link="logit"), quick = F, perm.parallel=T, nlambda=100, clus="id")
{
  
  cat("Requires libraries: ElasticLib, doParallel \n ") 
#  var.stand <- All_crossVars; 
  
#  mean.vec<-apply(data_l[,var.stand],2,mean)
#  sigma.vec<-apply(data_l[,var.stand],2,sd)
#  data_l[,var.stand]<-scale(data_l[,var.stand])
  
  # CHANGE
#  data_x  <- data_l[,c(clus,All_crossVars)] 
#  data_y <- data_l[,outcome]
  
  require(doParallel)
  
  rnd = list(id=~1)
  p1 <- length(clus)
  p <- ncol(data_x)-p1 # variables other than clustering indicators
  
  data_s <- cbind.data.frame(data_y, data_x) 
  CovarNames<- colnames(data_x)
  data <- WideData(data_s)
  
  T.max <- max(data$time)
  L.min <- min(data$time)
  time.length <- (T.max - L.min+1)

  t1 <- grep(paste0("X.",(T.max+1)),colnames(data))
  t2 <- grep(paste0("X.",1),colnames(data))
  
  
  colnames(data)[t2[1]:(t1[1]-1)] <- paste0("time.",1:time.length)
  colnames(data)[t1[1]:(t1[1]+(p-1))] <- colnames(data_x)[-1]
  
  covars <-  colnames(data_x)[-1]
  formula.x <- as.formula(paste("y~",paste0(covars,collapse="+"),sep=""))
  #formula.x <- as.formula(paste("y~",paste0(paste0("as.factor(",covars[1:2],")"),collapse="+"),"+",paste0(covars[-(1:2)],collapse="+"),sep=""))
  formula.time <- as.formula("~-1+time")
  if(is.null(index)){index <- 1:length(covars)}
#  cat("index =", index, "\n")
  
  
  n_g <- nrow(distinct(as_tibble(data[,clus,drop=F])))
  
  if(is.null(lambda.large)){warning("No lambda.large is specified.A default value 200 is assigned"); lambda.large=200}
  lambda <- lambda.large

  BIC_vec<-rep(Inf,length(lambda))
  
  nb =0 ; s=1;
  if(s > 0){nb <- n_g*s}
  
  if(is.null(init.theta)){Deltama.glm1<-as.matrix(t(rep(0,p+1+nb)))}else{Deltama.glm1<- init.theta}
  if(is.null(init.baseline)){Smooth.glm1<-as.matrix(t(rep(0,time.length)))}else{Smooth.glm1<- init.baseline}
   Q.glm1<-1e-1
   Q.course<-numeric()
   
  j<-1
  # Trying the estimate with large lambda 
  suppressWarnings(glm1 <- try(glmmElastic(formula.x,  
                          rnd = rnd,family = family, data = data, lambda=lambda[j], alpha=al,final.re=F,switch.NR=F,
                          control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                     diff.ord=0,penal= nus ,start=Smooth.glm1[j,]),
                                         method.final="EM", print.iter=F,print.iter.final=F, 
                                         eps.final=1e-4,epsilon=1e-4,complexity="non.zero",
                                         start=Deltama.glm1[j,],index=index, q_start=Q.glm1[j])),silent=T)) 
  if(class(glm1)[1]=="try-error")print("Red Flag 1 observed")
  lambda <- c(lambda,seq(from=glm1$lambda.max*1.1,to=0,length.out=nlambda))
  lambda <- lambda[-(nlambda+1)] 
  lambda.perm.vec  <- rep(NA, (lambda.perm.length+1))
  lambda.perm.vec[1] <- lambda.large
  Q.course.temp<-rep(0,length(lambda))
  
  #  Coef.ma.temp<-matrix(0,nrow=length(lambda),ncol=100)
  
  BIC_vec[j]<-glm1$bic
  # Coef.ma.temp[j,]<-glm1$coef[2:101]
  Deltama.glm1<-rbind(Deltama.glm1,glm1$Deltamatrix[glm1$conv.step,c(1:(p+1),(p+1+time.length+1):(p+1+time.length+nb))])
  Smooth.glm1<-rbind(Smooth.glm1,glm1$Deltamatrix[glm1$conv.step,c((p+1+1):(p+1+time.length))])
  Q.course.temp[j]<-glm1$StdDev
  Q.glm1<-c(Q.glm1,glm1$Q_long[[glm1$conv.step+1]])
  
  #=================================================  BIC loop =======================================================================================
  
  if(quick & Measure[1]!="perm"){stop("quick = T must be accompanied by \"perm\" measure. Option ignored"); quick = F}
  
  if(!quick){
   for(j in 2:length(lambda))
  {
    #print(paste("run",j,sep=": "))
     suppressWarnings( glm1 <- try(glmmElastic(formula.x,  
                            rnd = rnd,family = family, data = data, lambda=lambda[j],alpha=al,final.re=F,switch.NR=F,
                            control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                       diff.ord=0,penal=nus,start=Smooth.glm1[j,]),
                                           method.final="EM", print.iter=F,print.iter.final=F, 
                                           eps.final=1e-4,epsilon=1e-4,complexity="non.zero",
                                           start=Deltama.glm1[j,],index=index, q_start=Q.glm1[j])), silent=T)) 
  if(class(glm1)[1]=="try-error")
  {
    suppressWarnings(glm1 <- try(glmmElastic(formula.x,  
                            rnd = rnd,family = family, data = data, lambda=lambda[j],alpha=al,final.re=F,switch.NR=F,
                            control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                       diff.ord=0,penal=nus,start=Smooth.glm1[1,])
                                           ,method.final="EM", print.iter=T, print.iter.final=F,
                                           epsilon=1e-4,eps.final=1e-4,complexity="non.zero",
                                           start=Deltama.glm1[1,],index=index, q_start=Q.glm1[j])),silent=T)) 

    if(class(glm1)[1]!="try-error")
    {
      Deltama.glm1[j,]<-Deltama.glm1[1,]
      Q.glm1[j]<- Q.glm1[1]
      Smooth.glm1[j,]<-Smooth.glm1[1,]
    }  
  }
  
  
  if(class(glm1)[1]!="try-error")
  {
    BIC_vec[j]<-glm1$bic
    Deltama.glm1<-rbind(Deltama.glm1,glm1$Deltamatrix[glm1$conv.step,c(1:(p+1),(p+1+time.length+1):(p+1+time.length+nb))])
    Smooth.glm1<-rbind(Smooth.glm1,glm1$Deltamatrix[glm1$conv.step,c((p+1+1):(p+1+time.length))])
    Q.course.temp[j]<-glm1$StdDev
    Q.glm1<-c(Q.glm1,glm1$Q_long[[glm1$conv.step+1]])
      #  xx <- rbind(xx,c(glm1$aic,glm1$bic,glm1$loglik))
  }else{
    print(j)
    print("Red Flag 2 observed")
    Deltama.glm1<-rbind(Deltama.glm1,Deltama.glm1[j-1,])
    Smooth.glm1<-rbind(Smooth.glm1,Smooth.glm1[j-1,])
    Q.glm1<-c(Q.glm1,Q.glm1[j-1])
 #   xx <- rbind(xx,xx[j-1,])
  }
  
   }
    j_opt <- which.min(BIC_vec);
    
  }
#==========================================   Permutation method loop ========================================================
  perm.seed = 1234;
  if( Measure[1]=="perm"){
  if(perm.parallel){
      registerDoParallel(cores=future::availableCores());
    unpen <- c(1:ncol(data_y), ncol(data_y)+ (1:length(clus)), ncol(data_y)+length(clus)+which(is.na(index)));
    lambda.perm.vec.par <- foreach(j =  2:length(lambda.perm.vec), .combine=c) %dopar% {
      source("Preload.R"); sourceDir("./ElasticLib");
      mat = NA;
      set.seed(perm.seed+j);
      data.rnd <- data
      data.rnd[,unpen] <- data[sample.int(nrow(data)),unpen] 
      suppressWarnings(glm1.perm <- try(glmmElastic(formula.x,
                                   rnd = rnd,family = family, data = data.rnd, lambda=lambda.large, alpha=al,final.re=F,switch.NR=F,
                                   control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                              diff.ord=0,penal=nus,start=Smooth.glm1[1,]),
                                                  method.final="EM", print.iter=F,print.iter.final=F,
                                                  eps.final=1e-4,epsilon=1e-4,complexity="non.zero",
                                                  start=Deltama.glm1[1,],index=index,q_start=Q.glm1[1])),silent=T))
                                 if(class(glm1.perm)[1]!="try-error"){
                                 mat<-glm1.perm$lambda.max
                                }else{print("Red Flag 3 observed")}}
  #                         closeAllConnections()
                           lambda.perm.vec[-1] <- lambda.perm.vec.par
                           }else
                          {
                             unpen <- c(1:ncol(data_y), ncol(data_y)+ (1:length(clus)), ncol(data_y)+length(clus)+which(is.na(index))); 
                          for(j in 2:length(lambda.perm.vec)){
                                  set.seed(perm.seed+j);
                                  data.rnd <- data
                                  data.rnd[,unpen] <- data[sample.int(nrow(data)),unpen] 
                                  suppressWarnings(glm1.perm <- try(glmmElastic(formula.x,  
                                  rnd = rnd,family = family, data = data.rnd, lambda=lambda.large, alpha=al,final.re=F,switch.NR=F,
                                   control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                            diff.ord=0,penal=nus,start=Smooth.glm1[1,]),
                                                  method.final="Bres", print.iter=F,print.iter.final=F, 
                                                eps.final=1e-4,epsilon=1e-4,complexity="non.zero",
                                                start=Deltama.glm1[1,],index=index, q_start=Q.glm1[1])),silent=T))
    
                                   if(class(glm1.perm)[1]!="try-error"){
                                                                  lambda.perm.vec[j]<-glm1.perm$lambda.max
                                                                    }}}}
#===============================================================================================================
if(Measure[1]== "perm"){if(quick){j_opt <- 1}; lambda_opt <- median(lambda.perm.vec[-1], na.rm=T)}else{ lambda_opt=lambda[j_opt]} 
                      
  suppressWarnings(glm1.final <- try(glmmElastic(formula.x,  
                        rnd = rnd,family = family, data = data, lambda=lambda_opt,alpha=al,final.re=T,switch.NR=F,
                        control = list(smooth=list(formula=formula.time,nbasis=time.length,spline.degree=1,
                                                   diff.ord=0,penal=nus,start=Smooth.glm1[j_opt,]),
                                       method.final="EM", print.iter=T,print.iter.final=F, 
                                       eps.final=1e-4,epsilon=1e-4,complexity="non.zero",
                                       start=Deltama.glm1[j_opt,],index=index, q_start=Q.glm1[j_opt])),silent=T))
if(class(glm1.final)[1]=="try-error")print("Red Flag 4 observed")
cat("\n")  
xx <- data.frame(aic= glm1.final$aic, bic=glm1.final$bic, logLik=glm1.final$loglik) ;

sum <- cbind.data.frame(Method = Measure[1], quick =   quick, lambda= lambda_opt, alpha=al, nus=nus,xx) 
theta <- rbind(Deltama.glm1[j_opt,], c(glm1.final$coefficients, glm1$ranef))
time <- rbind(Smooth.glm1[j_opt,], glm1.final$smooth)
Q <- rbind(Q.glm1[j_opt], glm1.final$StdDev)
#=======================================================output to be returned ====================
out<- list(sum=sum,theta=theta, time=time, Q=Q, coefficients=glm1.final$coefficients, StdError=glm1.final$fixerror, baseline=glm1.final$smooth)
class(out) <- "Discnet"
return(out)
#================================================================================================
}



    