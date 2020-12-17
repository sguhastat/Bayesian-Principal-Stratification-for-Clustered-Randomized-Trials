rm(list=ls())
library(mvtnorm)
library(truncnorm)

## Get Data
load(file = "DataFile.Rdata")
y=datlist$y       ## outcome
X=datlist$X       ## predictor matrix
X <- as.matrix(X) 
N=datlist$N       ## sample size
M=datlist$M       ## no. of strata
treatment=datlist$treatment ## treatment status
D=datlist$D       ## intermediate variable
p.1=datlist$p.1   ## no. of cluster specific predictors
p.2=datlist$p.2   ## no. of individual specific predictors
K=datlist$K       ## no. of clusters
cluster.membership = datlist$cluster.membership
clust.index <- datlist$cluster  ## cluster membership indices
M <- 3  ## no. of strata
### Model fitting ##
## Three strata: complier = 0; always taker = 1; never taker = 2 ##

ind.clust <- list()
for(k in 1:K){ind.clust[[k]] <- which(clust.index==k)}         ## sample indices belonging to cluister k
z.1 <- sapply(1:N,function(ll){ifelse(treatment[ll]==1,1,0)})  ## treatment status
clust.tr <- ind.clust

## observed classes: (treatment or z.1, W), W is the intervention variable
W <- D
class.1 <- intersect(which(z.1==0),which(W==0))  ## compliers & never takers
class.2 <- intersect(which(z.1==0),which(W==1))  ## always takers
class.3 <- intersect(which(z.1==1),which(W==0))  ## never takers
class.4 <- intersect(which(z.1==1),which(W==1))  ## compliers & always takers

## class.1: 0+2, class.2: 1, class.3: 2, class.4: 0+1
#####################################################

## parameters
niter <- 50000           ## number of MCMC iterations
alpha_S.store <- list()  ## fixed effects in the data model for different strata 
beta_S.store <- list()   ## fixed effects in the weight model for different strata
for(m in 1:(M-1)){
   alpha_S.store[[m]] <- matrix(NA,niter,p.1+p.2)
   beta_S.store[[m]]  <- matrix(NA,niter,p.1+p.2)
}
alpha_S.store[[M]] <- matrix(NA,niter,p.1+p.2)
alpha_0.store <- matrix(NA,niter,M)  ## random intercept in the outcome model corresponding to treatment=0 for different strata
alpha_1.store <- matrix(NA,niter,M)  ## random intercept in the outcome model corresponding to treatment=1 for different strata
beta.1.store  <- matrix(NA,niter,K)  ## random intercept in the weight model 
beta.2.store  <- matrix(NA,niter,K)  ## random intercept in the weight model
S.store       <- matrix(NA,niter,N)  ## strata membership index for different MCMC iterates
omega.store   <- matrix(NA,niter,N)  ## latent variables corresponding to each sample to carry out Gibbs sampler for probit regression
G.1.store     <- matrix(NA,niter,N)  ## latent variables corresponding to each sample for weights
G.2.store     <- matrix(NA,niter,N)  ## latent variables corresponding to each sample for weights

## initialize
alpha_S <- list()
beta_S <- list()
for(m in 1:(M-1)){
   alpha_S[[m]] <- rnorm(p.1+p.2,0,0.1)
   beta_S[[m]]  <- rnorm(p.1+p.2,0,0.1)
}
alpha_S[[M]] <- rnorm(p.1+p.2,0,0.1)
alpha_0 <- rnorm(M)
alpha_1 <- rnorm(M)
beta.1  <- rnorm(K)
beta.2  <- rnorm(K)
S       <- sample(0:(M-1),N,replace=T,prob=rep(1/M,M))
S[class.1] <- sapply(1:length(class.1),function(ll){sample(c(0,2),1,prob=rep(1/2,2))}) 
S[class.2] <- rep(1,length(class.2))    ## the strata membership index is fixed for class.2
S[class.3] <- rep(2,length(class.3))    ## the strata membership index in fixed for class.3
S[class.4] <- sapply(1:length(class.4),function(ll){sample(c(0,1),1,prob=rep(1/2,2))})

omega   <- numeric()
G.1 <- rnorm(N,beta.1[clust.index]+c(X%*%beta_S[[1]]),1)
G.2 <- rnorm(N,beta.2[clust.index]+c(X%*%beta_S[[1]]),1)

## MCMC 

for(it in 1:niter){
   ## update omega (latent variable corresponding to each sample)
   index.1 <- c()
   index.2 <- c()
   index.3 <- c()
   index.4 <- c()
   for(j in 1:N){
      mem.j <- which(c(0:(M-1))==S[j])
      if((y[j]==1)&&(z.1[j]==0)){
         mean.1   <- alpha_0[mem.j]+sum(c(alpha_S[[mem.j]])*c(X[j,]))
         omega[j] <- rtruncnorm(1,0,Inf,mean.1,1)
         index.1  <- c(index.1,j)
      }else if((y[j]==1)&&(z.1[j]==1)){
         mean.1   <- alpha_1[mem.j]+sum(c(alpha_S[[mem.j]])*c(X[j,]))
         omega[j] <- rtruncnorm(1,0,Inf,mean.1,1)
         index.2  <- c(index.2,j)
      }else if((y[j]==0)&&(z.1[j]==0)){
         mean.1   <- alpha_0[mem.j]+sum(c(alpha_S[[mem.j]])*c(X[j,]))
         omega[j] <- rtruncnorm(1,-Inf,0,mean.1,1)
         index.3  <- c(index.3,j)
      }else{
         mean.1   <- alpha_1[mem.j]+sum(c(alpha_S[[mem.j]])*c(X[j,]))
         omega[j] <- rtruncnorm(1,-Inf,0,mean.1,1)
         index.4  <- c(index.4,j)
      }
    }

    ## update alpha_S (fixed effect vector corresponding to each stratum for the outcome model)
    for(m in 1:M){
        index.S.m   <- which(S==(m-1))
        X.S.m       <- X[index.S.m,]
        omega.S.m   <- omega[index.S.m]
        index.trt.0 <- which(treatment[index.S.m]==0)
        index.trt.1 <- which(treatment[index.S.m]==1)
        inter.S.m   <- rep(NA,length(index.S.m))
        inter.S.m[index.trt.0] <- rep(alpha_0[m],length(index.trt.0))
        inter.S.m[index.trt.1] <- rep(alpha_1[m],length(index.trt.1))
        var.alpha_S  <- chol2inv(chol(t(X.S.m)%*%X.S.m+diag(p.1+p.2)))
        mean.alpha_S <- var.alpha_S%*%t(X.S.m)%*%c(omega.S.m-inter.S.m)
        alpha_S[[m]] <- rmvnorm(1,mean.alpha_S,var.alpha_S)
    }

    ## update alpha_0 (random effect corresponding to each stratum and corresponding to treatment = 0 for the outcome model)
    for(m in 1:M){
       ind.m.0 <- intersect(which(S==(m-1)),which(z.1==0))
       if(length(ind.m.0)>1){
         omega.m.0 <- omega[ind.m.0]
         X.m.0     <- X[ind.m.0,]
         S.m.0     <- S[ind.m.0]
         slp.S.m.0  <- sapply(1:length(S.m.0),function(mm){sum(alpha_S[[which(c(0:(M-1))==S.m.0[mm])]]*X.m.0[mm,])})
         var.clust  <- 1/(length(ind.m.0)+1)
         mean.clust <- sum(omega.m.0-slp.S.m.0)*var.clust
         alpha_0[m] <- rnorm(1,mean.clust,sqrt(var.clust))
       }else if(length(ind.m.0)==1){
         omega.m.0 <- omega[ind.m.0]
         X.m.0     <- X[ind.m.0,]
         S.m.0     <- S[ind.m.0]
         slp.S.m.0  <- sum(alpha_S[[which(c(0:(M-1))==S.m.0)]]*X.m.0)
         var.clust  <- 1/(length(ind.m.0)+1)
         mean.clust <- sum(omega.m.0-slp.S.m.0)*var.clust
         alpha_0[m] <- rnorm(1,mean.clust,sqrt(var.clust))
        }else{
         alpha_0[m] <- rnorm(1)
        }      
    }

    ## update alpha_1 (random effect corresponding to each stratum and corresponding to treatment = 1 for the outcome model)
    for(m in 1:M){
       ind.m.1 <- intersect(which(S==(m-1)),which(z.1==1))
       if(length(ind.m.1)>1){
         omega.m.1 <- omega[ind.m.1]
         X.m.1     <- X[ind.m.1,]
         S.m.1     <- S[ind.m.1]
         slp.S.m.1  <- sapply(1:length(S.m.1),function(mm){sum(alpha_S[[which(c(0:(M-1))==S.m.1[mm])]]*X.m.1[mm,])})
         var.clust  <- 1/(length(ind.m.1)+1)
         mean.clust <- sum(omega.m.1-slp.S.m.1)*var.clust
         alpha_1[m] <- rnorm(1,mean.clust,sqrt(var.clust))
       }else if(length(ind.m.1)==1){
         omega.m.1 <- omega[ind.m.1]
         X.m.1     <- X[ind.m.1,]
         S.m.1     <- S[ind.m.1]
         slp.S.m.1  <- sum(alpha_S[[which(c(0:(M-1))==S.m.1)]]*X.m.1)
         var.clust  <- 1/(length(ind.m.1)+1)
         mean.clust <- sum(omega.m.1-slp.S.m.1)*var.clust
         alpha_1[m] <- rnorm(1,mean.clust,sqrt(var.clust))
        }else{
         alpha_1[m] <- rnorm(1)
        }      
    }


    ## update beta_S (fixed effect vector corresponding to each stratum for the weight model)
    S.1.f <- G.1-beta.1[clust.index]
    S.2.f <- G.2-beta.2[clust.index]
    v.beta.S.1 <- chol2inv(chol(t(X)%*%X+diag(p.1+p.2)))
    m.beta.S.1 <- v.beta.S.1%*%t(X)%*%c(S.1.f)
    m.beta.S.2 <- v.beta.S.1%*%t(X)%*%c(S.2.f)
    beta_S[[1]] <- rmvnorm(1,m.beta.S.1,v.beta.S.1)
    beta_S[[2]] <- rmvnorm(1,m.beta.S.2,v.beta.S.1)

    ## update beta.1 and beta.2 
    for(k in 1:K){
       resp.k.1 <- G.1[clust.tr[[k]]]-X[clust.tr[[k]],]%*%c(beta_S[[1]])
       resp.k.2 <- G.2[clust.tr[[k]]]-X[clust.tr[[k]],]%*%c(beta_S[[2]])
       beta.1[k] <- rnorm(1,sum(resp.k.1)/(length(clust.tr[[k]])+1),sqrt(1/(length(clust.tr[[k]])+1)))
       beta.2[k] <- rnorm(1,sum(resp.k.2)/(length(clust.tr[[k]])+1),sqrt(1/(length(clust.tr[[k]])+1)))
    }
    
    ## update G.1 and G.2
    for(j in 1:N){
       if(S[j]==0){
          mean.1 <- beta.1[clust.index[j]]+sum(X[j,]*c(beta_S[[1]]))
          mean.2 <- beta.2[clust.index[j]]+sum(X[j,]*c(beta_S[[2]]))
          G.1[j] <- rtruncnorm(1,-Inf,0,mean.1,1)
          G.2[j] <- rnorm(1,mean.2,1)
       }else if(S[j]==1){
          mean.1 <- beta.1[clust.index[j]]+sum(X[j,]*c(beta_S[[1]]))
          mean.2 <- beta.2[clust.index[j]]+sum(X[j,]*c(beta_S[[2]]))
          G.1[j] <- rtruncnorm(1,0,Inf,mean.1,1)
          G.2[j] <- rtruncnorm(1,-Inf,0,mean.2,1)
       }else{
          mean.1 <- beta.1[clust.index[j]]+sum(X[j,]*c(beta_S[[1]]))
          mean.2 <- beta.2[clust.index[j]]+sum(X[j,]*c(beta_S[[2]]))
          G.1[j] <- rtruncnorm(1,0,Inf,mean.1,1)
          G.2[j] <- rtruncnorm(1,0,Inf,mean.2,1)
       }
     }


    ## update S (strata membership index for each sample)
    for(j in 1:N){
       term1.1 <- log(pnorm(sapply(1:M,function(mm){alpha_0[mm]+sum(X[j,]*alpha_S[[mm]])})))
       term1.2 <- log(pnorm(sapply(1:M,function(mm){alpha_1[mm]+sum(X[j,]*alpha_S[[mm]])})))
       term2.1 <- beta.1[clust.index[j]]+sum(X[j,]*c(beta_S[[1]]))
       term2.2 <- beta.2[clust.index[j]]+sum(X[j,]*c(beta_S[[2]]))
       if(length(intersect(j,class.1))>0){     
          prob.1 <- log(pnorm(-term2.1))
          prob.2 <- log(pnorm(term2.1))+log(pnorm(term2.2))
          term2  <- c(prob.1,prob.2)
          ## for class 1, strata membership can be 0 (complier) or 2 (never taker)
          S[j]   <- ifelse(z.1[j]==0,sample(c(0,2),1,prob=exp(term1.1[c(1,3)]+term2)),sample(c(0,2),1,prob=exp(term1.2[c(1,3)]+term2))) 
       }else if(length(intersect(j,class.4))>0){        
          prob.1 <- log(pnorm(-term2.1))
          prob.2 <- log(pnorm(term2.1))+log(pnorm(-term2.2))
          term2  <- c(prob.1,prob.2)
          ## for class 4, strata membership can be 0 (complier) or 1 (always taker)
          S[j]   <- ifelse(z.1[j]==0,sample(c(0,1),1,prob=exp(term1.1[c(1,2)]+term2)),sample(c(0,1),1,prob=exp(term1.2[c(1,2)]+term2)))
       }else if(length(intersect(j,class.2))>0){
          S[j] <- 1  ## for class 2, strata membership is fixed (always taker)
       }else{
          S[j] <- 2  ## for class 3, strata membership is fixed (never taker)
       }
    }
    
    ## Store parameters
    for(m in 1:M){
       alpha_S.store[[m]][it,] <- alpha_S[[m]]
    }
    for(m in 1:(M-1)){
       beta_S.store[[m]][it,]  <- beta_S[[m]]
    }
    alpha_0.store[it,] <- alpha_0
    alpha_1.store[it,] <- alpha_1
    beta.1.store[it,]  <- beta.1
    beta.2.store[it,]  <- beta.2
    G.1.store[it,]     <- G.1
    G.2.store[it,]     <- G.2
    omega.store[it,]   <- omega
    S.store[it,]       <- S
    print(c(it,alpha_0))
}

####################
### Post-Processing ###
####################

burnin <- 40000
colMeans(alpha_S.store[[1]][(burnin+1):niter,])
colMeans(alpha_S.store[[2]][(burnin+1):niter,])
colMeans(alpha_S.store[[3]][(burnin+1):niter,])
colMeans(beta_S.store[[1]][(burnin+1):niter,])
colMeans(beta_S.store[[2]][(postburn+1):niter,])
colMeans(beta.1.store[(burnin+1):niter,])

### treatment effect calculation

prob.S.out <- array(NA,dim=c(niter-burnin,N,M))
treat.eff <- array(NA,dim=c(niter-burnin,N,M))
for(i in 1:(niter-burnin)){
  for(j in 1:N){
    term2.1 <- beta.1.store[i+burnin,clust.index[j]]+sum(X[j,]*c(beta_S.store[[1]][i+burnin,]))
    term2.2 <- beta.2.store[i+burnin,clust.index[j]]+sum(X[j,]*c(beta_S.store[[2]][i+burnin,]))
    prob.1 <- log(pnorm(-term2.1))
    prob.2 <- log(pnorm(term2.1))+log(pnorm(-term2.2))
    prob.3 <- log(pnorm(term2.1))+log(pnorm(term2.2))
    prob.S.out[i,j,]  <- exp(c(prob.1,prob.2,prob.3))
  }
  for(m in 1:M){
    treat.eff[i,,m] <- pnorm(alpha_1.store[i+burnin,m]+X%*%c(alpha_S.store[[m]][i+burnin,]))-pnorm(alpha_0.store[i+burnin,m]+X%*%c(alpha_S.store[[m]][i+burnin,]))
  }
}

treat.new <- matrix(NA,niter-burnin,N)
for(j in 1:N){
  if(length(intersect(j,class.1))>0){
    treat.new[,j] <- rowSums(treat.eff[,j,c(1,3)]*prob.S.out[,j,c(1,3)])/rowSums(prob.S.out[,j,c(1,3)])
  }else if(length(intersect(j,class.4))>0){
    treat.new[,j] <- rowSums(treat.eff[,j,c(1,2)]*prob.S.out[,j,c(1,2)])/rowSums(prob.S.out[,j,c(1,2)])
  }else if(length(intersect(j,class.2))>0){
    treat.new[,j] <- treat.eff[,j,2]
  }else{
    treat.new[,j] <- treat.eff[,j,3]
  }
}
overall.eff <- rowMeans(treat.new)
plot(density(overall.eff[seq(1,niter-burnin,by=4)]))
quantile(overall.eff[seq(1,niter-burnin,by=4)],c(0.025,.5,.975))

## Check convergence using trace plots of randomly selected parameters

plot(alpha_S.store[[1]][(burnin+1):niter,1], type = "l")
plot(alpha_S.store[[1]][(burnin+1):niter,4], type = "l")
plot(alpha_S.store[[2]][(burnin+1):niter,1], type = "l")
plot(alpha_S.store[[3]][(burnin+1):niter,4], type = "l")
dim(treat.new)

### Produce Final Plots
library(ggplot2)

## Overall Effect
df <- data.frame(TreatmentEffect=overall.eff[seq(1,niter-burnin,by=4)])
lower25 <- quantile(overall.eff[seq(1,niter-burnin,by=4)], probs = 0.025)  ## 2.5\% quantile of treatment effect
upper975 <- quantile(overall.eff[seq(1,niter-burnin,by=4)], probs = 0.975) ## 97.5% quantile of the treatment effect

head(df)
ggplot(df, aes(x=TreatmentEffect)) + geom_density(alpha=0.3,color="darkblue",fill="lightblue")+
  geom_vline(xintercept = lower25, col = 'coral', size = 1, linetype = "dashed")+
  geom_vline(xintercept = upper975, col = 'coral', size = 1, linetype = "dashed")+
  xlab("Overall Treatment Effect") + ylab("Density")+ ggtitle("Posterior Distribution of Treatment Effect")
theme(axis.text=element_text(size = 40), axis.title = element_text(size = 40, face = "bold"),axis.text.x = element_text(size = 40))

## ITT_W
ITT_W <- numeric()
for(i in 1:(niter-burnin)){
  ITT_W[i] <- length(which(S.store[i+burnin,]==0))/N
}
quantile(ITT_W, probs = c(0.025, 0.975))
mean(ITT_W)

## Proportion of Always-Takers
AT <- c()
for(i in 1:(niter-burnin)){
  AT[i] <- length(which(S.store[i+burnin,]==1))/N
}
quantile(AT, probs = c(0.025, 0.975))
mean(AT)

## Proportion of Never-Takers
NT <- c()
for(i in 1:(niter-burnin)){
  NT[i] <- length(which(S.store[i+burnin,]==2))/N
}
quantile(NT, probs = c(0.025, 0.975))
mean(NT)

strat.0 <- lapply(1:(niter-burnin),function(ll){which(S.store[ll+burnin,]==0)})
strat.1 <- lapply(1:(niter-burnin),function(ll){which(S.store[ll+burnin,]==1)})
strat.2 <- lapply(1:(niter-burnin),function(ll){which(S.store[ll+burnin,]==2)})

treat.complier <- sapply(1:(niter-burnin),function(ll){mean(treat.eff[ll,strat.0[[ll]],1])})
treat.always.taker <- sapply(1:(niter-burnin),function(ll){mean(treat.eff[ll,strat.1[[ll]],2])})
treat.never.taker <- sapply(1:(niter-burnin),function(ll){mean(treat.eff[ll,strat.2[[ll]],3])})
quantile(treat.complier,c(0.025,0.5,0.975))
quantile(treat.always.taker,c(0.025,0.5,0.975))
quantile(treat.never.taker,c(0.025,0.5,0.975))

## Complier Treatment Effect
dfc <- data.frame(TreatmentEffect=treat.complier)
lower25c <- quantile(treat.complier, probs = 0.025)
upper975c <- quantile(treat.complier, probs = 0.975)

head(dfc)
ggplot(dfc, aes(x=TreatmentEffect)) + geom_density(alpha=0.3,color="darkblue",fill="lightblue")+
  geom_vline(xintercept = lower25c, col = 'coral', size = 1, linetype = "dashed")+
  geom_vline(xintercept = upper975c, col = 'coral', size = 1, linetype = "dashed")+
  xlab("Complier Treatment Effect") + ylab("Density")+ ggtitle("Posterior Distribution of Treatment Effect")
theme(axis.text=element_text(size = 40), axis.title = element_text(size = 40, face = "bold"),axis.text.x = element_text(size = 40))

## Produce Plots for Always Taker and Never Taker - Treatment Effects
## Always Taker
dfa <- data.frame(TreatmentEffect=treat.always.taker)
lower25a <- quantile(treat.always.taker, probs = 0.025)
upper975a <- quantile(treat.always.taker, probs = 0.975)
ggplot(dfa, aes(x=TreatmentEffect)) + geom_density(alpha=0.3,color="darkblue",fill="lightblue")+
  geom_vline(xintercept = lower25a, col = 'coral', size = 1, linetype = "dashed")+
  geom_vline(xintercept = upper975a, col = 'coral', size = 1, linetype = "dashed")+
  xlab("Always Taker Treatment Effect") + ylab("Density")+ ggtitle("Posterior Distribution of Treatment Effect")
theme(axis.text=element_text(size = 40), axis.title = element_text(size = 40, face = "bold"),axis.text.x = element_text(size = 40))


dfn <- data.frame(TreatmentEffect=treat.never.taker)
lower25n <- quantile(treat.never.taker, probs = 0.025)
upper975n <- quantile(treat.never.taker, probs = 0.975)
ggplot(dfn, aes(x=TreatmentEffect)) + geom_density(alpha=0.3,color="darkblue",fill="lightblue")+
  geom_vline(xintercept = lower25n, col = 'coral', size = 1, linetype = "dashed")+
  geom_vline(xintercept = upper975n, col = 'coral', size = 1, linetype = "dashed")+
  xlab("Never Taker Treatment Effect") + ylab("Density")+ ggtitle("Posterior Distribution of Treatment Effect")
theme(axis.text=element_text(size = 40), axis.title = element_text(size = 40, face = "bold"),axis.text.x = element_text(size = 40))

## Posterior Means of Strata specific and Overall

meanall <- mean(overall.eff[seq(1,niter-burnin,by=4)])
meanc <- mean(treat.complier[seq(1,niter-burnin,by=4)])
meana <- mean(treat.always.taker[seq(1,niter-burnin,by=4)])
meann <- mean(treat.never.taker[seq(1,niter-burnin,by=4)])

meanall 
meanc 
meana 
meann 

## Finding the Coefficients - Outcome Model - alpha's
complier_coef <- apply(alpha_S.store[[1]][seq((burnin+1),niter,by=4),],2,quantile,c(0.025,.5,.975))
always_coef <- apply(alpha_S.store[[2]][seq((burnin+1),niter,by=4),],2,quantile,c(0.025,.5,.975))
never_coef  <- apply(alpha_S.store[[3]][seq((burnin+1),niter,by=4),],2,quantile,c(0.025,.5,.975))

complier_coef
complier_mean <- apply(alpha_S.store[[1]][seq((burnin+1),niter,by=4),],2,mean)
complier_mean
always_mean <- apply(alpha_S.store[[2]][seq((burnin+1),niter,by=4),],2,mean)
always_mean
always_coef
never_mean  <- apply(alpha_S.store[[3]][seq((burnin+1),niter,by=4),],2,mean)
never_mean
never_coef

## Finding the Coefficients - Model for inclusion probablities into a strata - beta's
beta1_coef <- apply(beta_S.store[[1]][seq((burnin+1),niter,by=4),],2,quantile,c(0.025,.5,.975))
beta1_coef
beta1_mean <- apply(beta_S.store[[1]][seq((burnin+1),niter,by=4),],2,mean)
beta1_mean

beta2_coef <- apply(beta_S.store[[2]][seq((burnin+1),niter,by=4),],2,quantile,c(0.025,.5,.975))
beta2_coef
beta2_mean <- apply(beta_S.store[[2]][seq((burnin+1),niter,by=4),],2,mean)
beta2_mean


