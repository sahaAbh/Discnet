#setwd("./<Curre directory >)
#setwd("D:/Codes")
#== The "Discnet.R" version  given below assumes one random effects due to indididual with gaussian distribution with vaeiance component , sigma^2. It can be easily modified to multiple effects with multiple variance components
#== Preload.R has auxilary functions to support the demo.R
#
#==Loading key R function, Discnet and all necessary libraries 
source("Preload.R")  # Required for demo.R
sourceDir("./ElasticLib") # Required for Discnet.R
source("Discnet.R") 

#==================================TRUE PARAMETERS FOR DATA GENERATION======
rho.vec= c(0.7,0.4) # Two types of pairwise correlations used
block=c(3,2) # Two blocks of compound symmetry used with correlations mentioned  
Lt_dist <-c(0.6,0.2,0.2) # Non-informative left truncation distribution,1 indicates no truncation
n <- 150    # sample size
p <- 150    # Number of covariates
cp <- 0.017 # Censor probability parameter that leads to roughly 20% censoring in this setting 
maxf <- 45 # maximum frequency  allowed for time-to-event  

baseline<-dgamma((1:10)-2, shape=5, scale = 1)*2- c(9,7,5,3,0,-1,-2,-4,-6,-8);#  Increasing baseline 
                                                                              #parameters with time
NZ= c(-4,-4,-4,8,8) #Vector of TRUE NON-ZERO covariate coefficients -
                    #starting from covariate 1 through the last NON-ZERO one
frailty_Q = 1 # frailty parameter: In this standard deviation of gaussian noise
T.max = 10 # number of baseline parameters
#=================== Generating data from a discrete frailty model==================================
seed= 124;
set.seed(seed)
Simout <- GenData(n=n,p=p,minfreq=2, maxfreq = maxf,baseline=baseline,cencor.prob=cp,rho.vec= rho.vec
                  ,block=block,frailty_Q=frailty_Q,Left.trun.prob=Lt_dist,NZ= NZ, T.max=T.max)


family= Simout$family
baseline=Simout$baseline
TrueTheta=Simout$beta 
test.seq=Simout$test.seq

data_s <- Simout$data_s
Censor <- 1-mean(data_s$delta)
cat("Observed Censoring:  ", Censor, "\n")
Trunc  <- mean(data_s$Lt >1) 
cat("Observed left truncation:  ", Trunc, "\n")

T.max <- max(data_s$time)
L.min <- min(data_s$time)
time.length <- (T.max - L.min+1)

outcome <- c("time", "delta", "Lt");
clus <- c("id") # id specifying subjects to include subject specific frailty
data_x  <- data_s[,c(clus,paste0("X.",1:p))] 
data_y <- data_s[,outcome]

#=================================  TUNING PEBNALTY PARAMETERS================
al.vec <- c(1,0.95,0.8,0.7,0.6,0.5) # grid for alpha (elastic net penalty) 
nus.vec=c(15,25,50,100) # grid for nu, penalty parameters for baselines
ParamTab <- CrossParamTab(c( "alpha","nus"), list(alpha=al.vec,nus=nus.vec)) 
                       # lists two way grids for alpha and nu
Mes <- "perm" # method to choose lambda (elastic net penalty),  
       # "bic" is the other alternative   
lambda.large= 200  # An initial choice of  large lasso penalty (of elastic net) that forces 
                   # all selection parameters to become 0 

tic1 <- Sys.time()
index <- 1:p; # NA indicates non-penalization, same group variables have same index 
# for example, from among V1,..V10 if only V1,V2 are not penalized; V3,V4,v5 are dummy codes 
# for a categorical variable with 4 levels whereas rest are conts then index=c(NA,NA,3,3,3,3,4,5,6,7)  
for(rowind in 1:nrow(ParamTab)){

al = ParamTab[rowind,"alpha"];
nus = ParamTab[rowind,"nus"]

print(paste0("seed=",seed,", al=",al,", nus=",nus))

Dspath <- Discnet(data_y, data_x, al=al,nus=nus, index=index,  Measure = Mes,  lambda.large= lambda.large)

if(rowind==1){Sum.tab <- cbind.data.frame(data.seed = seed, Dspath$sum)}else
{
  Sum.tab <- rbind.data.frame(Sum.tab, cbind.data.frame(data.seed = seed, Dspath$sum))  
  }
}
print(Sys.time()-tic1)
print(Sum.tab)

tune_opt <- which.min(Sum.tab$bic)
Dspath <- Discnet(data_y, data_x, al=Sum.tab$alpha[tune_opt],nus=Sum.tab$nus[tune_opt],index=index,  Measure = Mes, lambda.large= lambda.large)

theta_est <- Dspath$theta[2,2:(p+1),drop=F] # First coeffiecient is the intercept
theta_T <- TrueTheta[11:(11+p-1)]
theta_err <- sum((theta_T-theta_est)^2)

#---------False Posities and False Negatives
epsilon <- 0
NZ =  1*( abs(theta_T) > epsilon)[1:p]
Pmat <- 1*(abs(theta_est) > epsilon) 
FN <- sum(NZ)-Pmat %*% NZ
FP <-  Pmat %*% (!NZ)

Data.sum <-  cbind.data.frame(data.seed = seed,Censor=Censor,Trunc=Trunc, frailty=Dspath$Q[2],Dspath$sum, FN=FN,FP=FP,theta_err=theta_err)  #Spline case
print(Data.sum)
#======================= Summarizing the output=====
names(summary(Dspath))
summary(Dspath)$coefficients[1:10,] # NA means non-selection of the variable
summary(Dspath)$baseline.eff  #Baseline estimates
summary(Dspath)$Q             #frailty estimates

