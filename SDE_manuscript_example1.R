##############################################################################
##############################################################################
#This code is to conduct analyses for the example 1 in the manuscript titled
#Adding noise to Markov cohort models
#under review at MDM
#author: Rowan Iskandar
#email: rowan_iskandar@brown.edu
#last edited: May 10 2019
#The code for the numerical exercise is available under a GNU GPL license and can be found at https://github.com/rowaniskandar/SDE.
##############################################################################
##############################################################################
#loading relevant libraries
library(yuima)
library(deSolve)
library(ggplot2)
library(reshape2)
library(expm)
library(pracma)
library(markovchain)
library(matrixcalc)
library(openxlsx)
currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
set.seed(12345)
##############################################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################################################################
##############################################################################
#parameter values for example 1
#reference: Iskandar, Rowan. "A theoretical foundation for state-transition cohort models in health decision analysis." PloS one 13.12 (2018): e0205543.
##############################################################################
#parameter values for transition rates in the 4-state model (note: these are rates, not probabilities)
c12 <- 0.05
c13 <- 0.01
c14 <- 0.001
c23 <- 0.1
c24 <- 0.05
c34 <- 1

popsize <- 1000
ns <- 4 #number of states
nMC <-10000 # number of micro sim iterations
t0 <- 0 #initial time
tfin <-50 #final time
dt <- 1 #time step
tpoints <- seq(t0,tfin,dt) #number of discrete time points

params <- c(c12 ,
               c13 ,
               c14 ,
               c23 ,
               c24 ,
               c34 )
init <-c(popsize,0,0,0) #initial state configuration
#generator matrix
Q <-matrix(c(-(c12+c13+c14),0,0,0,c12,-(c23+c24),0,0,c13, c23, -(c34),0,c14, c24,c34,0), nrow=4,ncol=4)
tstep <- dt
##############################################################################
##############################################################################
##SDE 
####2019 03 25
#vector v
# prop1 <- function(n1){ return(c12*n1)}
# prop2 <- function(n1){ return(c13*n1)}
# prop3 <- function(n1){ return(c14*n1)}
# prop4 <- function(n2){ return(c23*n2)}
# prop5 <- function(n2){ return(c24*n2)}
# prop6 <- function(n3){ return(c34*n3)}
c <-c(c12,c13,c14,c23,c24,c34)
vprop <- function(n,c){
    v <-matrix(rep(0,length(c)),nrow=length(c),ncol=1,byrow = TRUE)
    v[1,1] <- c[1]*n[1]
    v[2,1] <- c[2]*n[1]
    v[3,1] <- c[3]*n[1]
    v[4,1] <- c[4]*n[2]
    v[5,1] <- c[5]*n[2]
    v[6,1] <- c[6]*n[3]
  return(v)
}

#matrid d=change
d <- matrix(c(-1,1,0,0,-1,0,1,0,-1,0,0,1,0,-1,1,0,0,-1,0,1,0,0,-1,1),nrow=6,ncol=4,byrow=TRUE)

#change mean
change_mean <-function(n,c,d){
  A <- as.vector(t(vprop(n,c)))%*%d
  return(A)
}

#change variance
change_var <-function(n,c,d){
  B<-matrix(rep(0,length(n)*length(n)),nrow=length(n), ncol=length(n))
  num_trans <- length(c)
  num_states <- length(n)
  for (i in 1:num_states){
    for (j in 1:num_states){
      sum=0
      for (k in 1:num_trans){
        sum=sum+vprop(n,c)[k,1]*d[k,i]*d[k,j]
      }
      B[i,j] <-sum
    }
  }  
  return(B)
}
#####euler maruyama for one iteration
SDE.trace <- matrix(rep(0,length(tpoints)*ns),nrow=length(tpoints),ncol=ns)
SDE.trace[1,] <- init

#####yuima package method
sol <- c("n1","n2","n3","n4")
A <-c("-n1*(c12+c13+c14)","n1*c12-n2*(c23+c24)","n1*c13+n2*c23-n3*c34","n1*c14+n2*c24+n3*c34")
B <-matrix(c("n1*(c12+c13+c14)","-n1*c12","-n1*c13","-n1*c14","-n1*c12","n1*c12+n2*c23+n2*c24","-n2*c23","-n2*c24","-n1*c13","-n2*c23","n1*c13+n2*c23+n3*c34","-n3*c34","-n1*c14","-n2*c24","-n3*c34","n1*c14+n2*c24+n3*c34") ,nrow=4,ncol=4)
#A <-c("-n1*(0.05+0.01+0.001)","n1*0.05-n2*(0.1+0.05)","n1*0.01+n2*0.1-n3*1","n1*0.001+n2*0.05+n3*1")
# B <-matrix(c("n1*(0.05+0.01+0.001)","-n1*0.05","-n1*0.01","-n1*0.001","-n1*0.05","n1*0.05+n2*0.1+n2*0.05","-n2*0.1","-n2*0.05","-n1*0.01","-n2*0.1","n1*0.01+n2*0.1+n3*1","-n3*1","-n1*0.001","-n2*0.05","-n3*1","n1*0.001+n2*0.05+n3*1") ,nrow=4,ncol=4)
MCtrace_all <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
MCtrace_all2 <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
tic()
for(m in 1:nMC){
  mod <- setModel(drift=A,diffusion=B,state.variable=sol,solve.variable = sol)
  samp <- setSampling(Initial=t0,Terminal=tfin, n=length(tpoints)-1)
  yui <-setYuima(model=mod,sampling=samp)
  sde.out <- simulate(yui,xinit=init,true.parameter=params)
  sde<-get.zoo.data(sde.out)
  MCtrace_all[m, ,1]<-as.vector(sde$`Series 1`)
  MCtrace_all[m, ,2]<-as.vector(sde$`Series 2`)
  MCtrace_all[m, ,3]<-as.vector(sde$`Series 3`)
  MCtrace_all[m, ,4]<-as.vector(sde$`Series 4`)
}
toc()
#data.sde<-as.data.frame(sde)
data.MCtrace_all <- as.data.frame(MCtrace_all)
dum1 <-colSums(MCtrace_all)
meanMC <-dum1/nMC
varMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdposMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdnegMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)

for (m in 1:nMC){
  MCtrace_all2[m, ,1] <-  (MCtrace_all[m, ,1]-meanMC[,1])^2
  MCtrace_all2[m, ,2] <-  (MCtrace_all[m, ,2]-meanMC[,2])^2
  MCtrace_all2[m, ,3] <-  (MCtrace_all[m, ,3]-meanMC[,3])^2
  MCtrace_all2[m, ,4] <-  (MCtrace_all[m, ,4]-meanMC[,4])^2
}

dum2 <-colSums(MCtrace_all2)
varMC <-dum2/nMC
sdposMC <- meanMC

sdposMC[,1] <- meanMC[,1] + sqrt(varMC[,1])
sdnegMC[,1] <- meanMC[,1] - sqrt(varMC[,1])
sdposMC[,2] <- meanMC[,2] + sqrt(varMC[,2])
sdnegMC[,2] <- meanMC[,2] - sqrt(varMC[,2])
sdposMC[,3] <- meanMC[,3] + sqrt(varMC[,3])
sdnegMC[,3] <- meanMC[,3] - sqrt(varMC[,3])
sdposMC[,4] <- meanMC[,4] + sqrt(varMC[,4])
sdnegMC[,4] <- meanMC[,4] - sqrt(varMC[,4])

sdposMC_d <- sdposMC
sdnegMC_d <- sdnegMC
meanMC_d <- meanMC

#sdposMC_d[1,] <- c(popsize,0,0,0)
#sdnegMC_d[1,] <- c(popsize,0,0,0)
#meanMC_d[1,] <- c(popsize,0,0,0)

for (t in 2:length(tpoints)){
  sdposMC_d[t,] <- sdposMC[t-1,]
  sdnegMC_d[t,] <- sdnegMC[t-1,]
  meanMC_d[t,] <- meanMC[t-1,]
}

data.SDEyui.mean <- data.frame(meanMC)
data.SDEyui.var <- data.frame(varMC)
data.SDEyui.sdpos <- data.frame(sdposMC)
data.SDEyui.sdneg <- data.frame(sdnegMC)


##############################################################################
##############################################################################
##ordinary differential equation version of cohort model (Equation 4 in the reference)
odemod <- function(t,n,parameters){
  # c12 <- 0.05
  # c13 <- 0.01
  # c14 <- 0.001
  # c23 <- 0.1
  # c24 <- 0.05
  # c34 <- 1
  list(c(-n[1]*(c12+c13+c14),n[1]*c12-n[2]*(c23+c24),n[1]*c13+n[2]*c23-n[3]*c34,n[1]*c14+n[2]*c24+n[3]*c34))
}
ode.out <- ode(y = init, func = odemod,times = 0:50, parms = 1)
data.ode <-as.data.frame(ode.out)
# ##############################################################################
# ##############################################################################
# ##Cohort simulation (Equation 3), the transition probabilities are calculated based on the exact solution of Kolmogorov equation (Welton and Ades 2005)
p11 <- function(params,t) {
  s1 <-params[1]+params[2]+params[3]
  p<-exp(-s1*t)
return(p)}
p12 <- function(params,t) {
  s1 <-params[1]+params[2]+params[3]
  s2 <-params[4]+params[5]
  p<-params[1]*(exp(-s2*t)-exp(-s1*t))/(s1-s2)
  return(p)}
p13 <- function(params,t) {
  s1 <-params[1]+params[2]+params[3]
  s2 <-params[4]+params[5]
  s3 <-params[6]
  p1<-params[2]*(exp(-s3*t)-exp(-s1*t))/(s1-s3)
  p2<-(params[1]*params[4]*((s1-s2)*exp(-s3*t)-(s1-s3)*exp(-s2*t)+(s2-s3)*exp(-s1*t))) /((s1-s2)*(s1-s3)*(s2-s3))
  p<-p1+p2
  return(p)}
p14 <- function(params,t) {
  p1 <- p11(params,t)
  p2 <- p12(params,t)
  p3 <- p13(params,t)
  p <- 1-p1-p2-p3
  return(p)}

p21 <- 0
p22 <-function(params,t) {
  s2 <-params[4]+params[5]
  p<-exp(-s2*t)
  return(p)}

p23 <-function(params,t) {
  s2 <-params[4]+params[5]
  s3 <- params[6]
  p<-params[4]*(exp(-s3*t)-exp(-s2*t))/(s2-s3)
  return(p)}
p24 <- function(params,t) {
  p1 <- 0
  p2 <- p22(params,t)
  p3 <- p23(params,t)
  p <- 1-p1-p2-p3
  return(p)}

p31 <- 0
p32 <- 0

p33 <-function(params,t) {
  s3 <-params[6]
  p<-exp(-s3*t)
  return(p)}

p34 <- function(params,t) {
  p <- 1-p33(params,t)
  return(p)}

p41 <- 0
p42 <- 0
p43 <- 0
p44 <- 1

prob <-function(params,t){
  P <- matrix(c(p11(params,t) ,p21,p31,p41,p12(params,t) ,p22(params,t) ,p32,p42,p13(params,t) ,p23(params,t) ,p33(params,t) ,p43,p14(params,t) ,p24(params,t) ,p34(params,t) ,p44 ),nrow=4,ncol=4)
  return(P)
  }
CM.trace <- matrix(rep(0,length(tpoints)*4),nrow=length(tpoints),ncol=4)
p.temp <- 0
for(i in 1:length(tpoints)){
  # for(j in 1:i){
  #   if(j==1){
  #     p.temp <- prob(params,1)
  #   }
  #   else{
  #     p.temp <- p.temp%*%prob(params,j)
  #   }
  # }
  #CM.trace[i,] <- init%*%p.temp
  CM.trace[i,] <- init%*%prob(params,i)
}

CM.trace_d <- CM.trace
CM.trace_d[1,] <-c(popsize, 0 ,0 ,0)
for (t in 2:length(tpoints)){
  CM.trace_d[t,] <- CM.trace[t-1,]
}

data.CM <- as.data.frame(CM.trace_d)
##############################################################################
##############################################################################
#microsimulation of the 4-state transition model
sickStates <- c("healthy","mild","severe","death")
trprob <- expm(Q*tstep)

sickP <-matrix(data = trprob,  nrow = ns, ncol=ns, dimnames = list(sickStates,sickStates))
sickMC <- new("markovchain", states = sickStates , byrow=TRUE, transitionMatrix =sickP , name = "sick")

N <- popsize
tfinal <- length(tpoints)
initialState <- c(1,0,0,0)

sicktrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
sickbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
healthtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
healthbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
mildtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
mildbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
severetrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
severebool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
deathtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
deathbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
MCtrace <- matrix(rep(0,length(tpoints)*ns), nrow = length(tpoints), ncol = ns)
MCtrace_all <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
MCtrace_all2 <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))

tic()
#start outer outer loop
for (m in 1:nMC){
  for (i in 1:N){
    sicktrace[i,] <- rmarkovchain(n = tfinal, object = sickMC, t0 = "healthy")
  }
  
  for (i in 1:N){
    for (j in 1:tfinal){
      if (sicktrace[i,j]=="healthy"){healthbool[i,j]=1}
      else {healthbool[i,j]=0}
      if (sicktrace[i,j]=="mild"){mildbool[i,j]=1}
      else {mildbool[i,j]=0}
      if (sicktrace[i,j]=="severe"){severebool[i,j]=1}
      else {severebool[i,j]=0}
      if (sicktrace[i,j]=="death"){deathbool[i,j]=1}
      else {deathbool[i,j]=0}
    }
  }
  MCtrace_all[m, ,1] <- colSums(healthbool, dims=1)
  MCtrace_all[m, ,2] <- colSums(mildbool, dims=1)
  MCtrace_all[m, ,3] <- colSums(severebool, dims=1)
  MCtrace_all[m, ,4] <- colSums(deathbool, dims=1)
}
toc()
dum1 <-colSums(MCtrace_all)
meanMC <-dum1/nMC
varMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdposMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdnegMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)

for (m in 1:nMC){
  MCtrace_all2[m, ,1] <-  (MCtrace_all[m, ,1]-meanMC[,1])^2
  MCtrace_all2[m, ,2] <-  (MCtrace_all[m, ,2]-meanMC[,2])^2
  MCtrace_all2[m, ,3] <-  (MCtrace_all[m, ,3]-meanMC[,3])^2
  MCtrace_all2[m, ,4] <-  (MCtrace_all[m, ,4]-meanMC[,4])^2
}

dum2 <-colSums(MCtrace_all2)
varMC <-dum2/nMC
sdposMC <- meanMC

sdposMC[,1] <- meanMC[,1] + sqrt(varMC[,1])
sdnegMC[,1] <- meanMC[,1] - sqrt(varMC[,1])
sdposMC[,2] <- meanMC[,2] + sqrt(varMC[,2])
sdnegMC[,2] <- meanMC[,2] - sqrt(varMC[,2])
sdposMC[,3] <- meanMC[,3] + sqrt(varMC[,3])
sdnegMC[,3] <- meanMC[,3] - sqrt(varMC[,3])
sdposMC[,4] <- meanMC[,4] + sqrt(varMC[,4])
sdnegMC[,4] <- meanMC[,4] - sqrt(varMC[,4])

sdposMC_d <- sdposMC
sdnegMC_d <- sdnegMC
meanMC_d <- meanMC

sdposMC_d[1,] <- c(popsize,0,0,0)
sdnegMC_d[1,] <- c(popsize,0,0,0)
meanMC_d[1,] <- c(popsize,0,0,0)

for (t in 2:length(tpoints)){
  sdposMC_d[t,] <- sdposMC[t-1,]
  sdnegMC_d[t,] <- sdnegMC[t-1,]
  meanMC_d[t,] <- meanMC[t-1,]
}

data.MC.mean <- data.frame(meanMC_d)
data.MC.var <- data.frame(varMC)
data.MC.sdpos <- data.frame(sdposMC_d)
data.MC.sdneg <- data.frame(sdnegMC_d)

##############################################################################
##############################################################################
#master equation method  
#mean of master equation (Equation 8)
# pgf_mean <- function(t){
#     C <-expm(t*Q)
#     sum_C <- 0
#     for (j in 1:ns){
#       sum_C <- sum_C + C[1,j]
#     }
#     #mean_first_term <- n0*(sum_C)^(n0-1)
#     mean_first_term <- n0
#     pgf1 <- mean_first_term*C[1,1]
#     pgf2 <- mean_first_term*C[1,2]
#     pgf3 <- mean_first_term*C[1,3]
#     pgf4 <- mean_first_term*C[1,4]
#     #print(sum_C)
#   return(c(pgf1,pgf2,pgf3,pgf4))
# }
# #variance of master equation (Equation 9)
# pgf_var <-function(t){
#   mean <- pgf_mean(t)
#   C <-expm(tpoints[i]*Q)
#   D <- Q%*%C
#   sum_C <- 0
#   for (j in 1:ns){
#     sum_C <- sum_C + C[1,j]
#   }
#   
#   var1 <- mean[1]*((n0-1)*(1/sum_C)*C[1,1]+1-mean[1])
#   var2 <- mean[2]*((n0-1)*(1/sum_C)*C[1,2]+1-mean[2])
#   var3 <- mean[3]*((n0-1)*(1/sum_C)*C[1,3]+1-mean[3])
#   var4 <- mean[4]*((n0-1)*(1/sum_C)*C[1,4]+1-mean[4])
#   
#   return(c(var1,var2,var3,var4))
# }
# 
# fac <- 0
# flag.true <- FALSE
# pdf_count <- function(state.config,t){
#     C <-expm(t*Q)
#     pdf <- c(0,0,0,0)
#     for (i in 1:ns){
#       x <- state.config[i]
#       fac <- factorial(x)
#       pdf[i] <- C[1,i]^(x)/fac
#     }
#     pdf_out <-popsize*pdf[1]*pdf[2]*pdf[3]*pdf[4]
#     return(c(pdf_out))
# }
# 
# count.data <- seq(1,popsize)
# count.design.frame <- expand.grid(Factor1 = count.data, Factor2 = count.data, Factor3 = count.data, Factor4 = count.data)
# count.design <- data.matrix(count.design.frame,rownames.force = NA)
# data.pdf <- matrix(rep(0,popsize^4),nrow =popsize^4,ncol=1 )
# 
# sum.pdf <- 0
# for (i in 1:popsize^4){
#   data.pdf[i]=pdf_count(count.design[i,], 1)
#   sum.pdf <- sum.pdf + data.pdf[i]
# }
# 
# for (i in 1:popsize^4){
#   data.pdf[i] <- data.pdf[i]/sum.pdf
# }
# 
# n0 <- init[1]
# pgf.trace <- matrix(rep(0,ns*length(tpoints)),length(tpoints),ns)
# pgf.var.trace <- matrix(rep(0,ns*length(tpoints)),length(tpoints),ns)
# pgf.sdpos.trace <- matrix(rep(0,ns*length(tpoints)),length(tpoints),ns)
# pgf.sdneg.trace <- matrix(rep(0,ns*length(tpoints)),length(tpoints),ns)
# 
# for (i in 1:length(tpoints)){
#   mean <- pgf_mean(tpoints[i])
#   var <- pgf_var(tpoints[i])
#   
#   pgf.trace[i,1] <- mean[1]
#   pgf.trace[i,2] <- mean[2]
#   pgf.trace[i,3] <- mean[3]
#   pgf.trace[i,4] <- mean[4]
#   
#   pgf.var.trace[i,1] <- var[1]
#   pgf.var.trace[i,2] <- var[2]
#   pgf.var.trace[i,3] <- var[3]
#   pgf.var.trace[i,4] <- var[4]
#   
#   pgf.sdpos.trace[i,1] <- mean[1]+sqrt(var[1])
#   pgf.sdpos.trace[i,2] <- mean[2]+sqrt(var[2])
#   pgf.sdpos.trace[i,3] <- mean[3]+sqrt(var[3])
#   pgf.sdpos.trace[i,4] <- mean[4]+sqrt(var[4])
#   
#   pgf.sdneg.trace[i,1] <- mean[1]-sqrt(var[1])
#   pgf.sdneg.trace[i,2] <- mean[2]-sqrt(var[2])
#   pgf.sdneg.trace[i,3] <- mean[3]-sqrt(var[3])
#   pgf.sdneg.trace[i,4] <- mean[4]-sqrt(var[4])
# }
# data.pgf <- as.data.frame(pgf.trace)
# data.pgf.var <- as.data.frame(pgf.var.trace)
# data.pgf.sdpos <- as.data.frame(pgf.sdpos.trace)
# data.pgf.sdneg <- as.data.frame(pgf.sdneg.trace)

# for (i in 1:length(tpoints)){
#   C <-expm(tpoints[i]*Q)
#   sum_C <- 0
#   for (j in 1:ns){
#     sum_C <- sum_C + C[1,j]
#   }
#   mean_first_term <- n0*(sum_C)^(n0-1)
#   pgf.var.trace[i,1] <- mean_first_term*C[1,1]
#   pgf.var.trace[i,2] <- mean_first_term*C[1,2]
#   pgf.var.trace[i,3] <- mean_first_term*C[1,3]
#   pgf.var.trace[i,4] <- mean_first_term*C[1,4]
# }
##############################################################################
##############################################################################
#combining data
# n1 <- cbind(data.ode$time,data.ode$`1`,data.CM$V1,data.pgf$V1, data.pgf.sdpos$V1, data.pgf.sdneg$V1, data.MC$X1)
# data.n1 <- as.data.frame(n1)
# n2 <- cbind(data.ode$time,data.ode$`2`,data.CM$V2,data.pgf$V2, data.pgf.sdpos$V2, data.pgf.sdneg$V2, data.MC$X2)
# data.n2 <- as.data.frame(n2)
# n3 <- cbind(data.ode$time,data.ode$`3`,data.CM$V3,data.pgf$V3, data.pgf.sdpos$V3, data.pgf.sdneg$V3, data.MC$X3)
# data.n3 <- as.data.frame(n3)
# n4 <- cbind(data.ode$time,data.ode$`4`,data.CM$V4,data.pgf$V4, data.pgf.sdpos$V4, data.pgf.sdneg$V4, data.MC$X4)
# data.n4 <- as.data.frame(n4)
#variance comparison
# n1 <- cbind(data.ode$time,data.ode$`1`,data.CM$V1,data.SDE.mean$X1, data.SDE.sdpos$X1, data.SDE.sdneg$X1, data.MC.mean$X1, data.MC.sdpos$X1, data.MC.sdneg$X1)
# data.n1 <- as.data.frame(n1)
# n2 <- cbind(data.ode$time,data.ode$`2`,data.CM$V2,data.SDE.mean$X2, data.SDE.sdpos$X2, data.SDE.sdneg$X2, data.MC.mean$X2, data.MC.sdpos$X2, data.MC.sdneg$X2)
# data.n2 <- as.data.frame(n2)
# n3 <- cbind(data.ode$time,data.ode$`3`,data.CM$V3,data.SDE.mean$X3, data.SDE.sdpos$X3, data.SDE.sdneg$X3, data.MC.mean$X3, data.MC.sdpos$X3, data.MC.sdneg$X3)
# data.n3 <- as.data.frame(n3)
# n4 <- cbind(data.ode$time,data.ode$`4`,data.CM$V4,data.SDE.mean$X4, data.SDE.sdpos$X4, data.SDE.sdneg$X4, data.MC.mean$X4, data.MC.sdpos$X4, data.MC.sdneg$X4)
# data.n4 <- as.data.frame(n4)
#yuima
n1yui <- cbind(data.ode$time,data.ode$`1`,data.CM$V1,data.SDEyui.mean$X1, data.SDEyui.sdpos$X1, data.SDEyui.sdneg$X1, data.MC.mean$X1, data.MC.sdpos$X1, data.MC.sdneg$X1)
data.n1yui <- as.data.frame(n1yui)
n2yui <- cbind(data.ode$time,data.ode$`2`,data.CM$V2,data.SDEyui.mean$X2, data.SDEyui.sdpos$X2, data.SDEyui.sdneg$X2, data.MC.mean$X2, data.MC.sdpos$X2, data.MC.sdneg$X2)
data.n2yui <- as.data.frame(n2yui)
n3yui <- cbind(data.ode$time,data.ode$`3`,data.CM$V3,data.SDEyui.mean$X3, data.SDEyui.sdpos$X3, data.SDEyui.sdneg$X3, data.MC.mean$X3, data.MC.sdpos$X3, data.MC.sdneg$X3)
data.n3yui <- as.data.frame(n3yui)
n4yui <- cbind(data.ode$time,data.ode$`4`,data.CM$V4,data.SDEyui.mean$X4, data.SDEyui.sdpos$X4, data.SDEyui.sdneg$X4, data.MC.mean$X4, data.MC.sdpos$X4, data.MC.sdneg$X4)
data.n4yui <- as.data.frame(n4yui)
##############################################################################
##############################################################################
#write table tab delimited
write.table(data.n1,paste(currpath,"n1.txt" ,sep=""), sep="\t")
write.table(data.n2,paste(currpath,"n2.txt" ,sep=""), sep="\t")
write.table(data.n3,paste(currpath,"n3.txt" ,sep=""), sep="\t")
write.table(data.n4,paste(currpath,"n4.txt" ,sep=""), sep="\t")

write.xlsx(data.n1,paste(currpath,"n1.xlsx" ,sep=""))
write.xlsx(data.n2,paste(currpath,"n2.xlsx" ,sep=""))
write.xlsx(data.n3,paste(currpath,"n3.xlsx" ,sep=""))
write.xlsx(data.n4,paste(currpath,"n4.xlsx" ,sep=""))

write.xlsx(data.n1,paste(currpath,"n1.xlsx" ,sep=""))
write.xlsx(data.n2,paste(currpath,"n2.xlsx" ,sep=""))
write.xlsx(data.n3,paste(currpath,"n3.xlsx" ,sep=""))
write.xlsx(data.n4,paste(currpath,"n4.xlsx" ,sep=""))

##############################################################################
##############################################################################
#plotting
plot.n1<-ggplot(data.n1, aes(V1,y=value,color=variable)) +  geom_line(aes(y=V2,color="ODE"))+ geom_line(aes(y=V3,color="CM"))+geom_line(aes(y=V4,color="SDE"))+geom_line(aes(y=V5,color="SDE+sd")) +geom_line(aes(y=V6,color="SDE-sd")) +geom_line(aes(y=V7,color="MSM"))+geom_line(aes(y=V8,color="MSM+sd"))+geom_line(aes(y=V9,color="MSM-sd"))+scale_shape_discrete(solid=T)  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x="Time",y=expression("Number of individuals in "*"S"[1])) 
  
ggsave(paste(currpath,"/n1.tiff",sep="") , units="in", width=6, height=5, dpi=1200, compression = 'lzw')

plot.n2<-ggplot(data.n2, aes(V1,y=value,color=variable)) +  geom_line(aes(y=V2,color="ODE")) + geom_line(aes(y=V3,color="CM"))+geom_line(aes(y=V4,color="SDE")) +geom_line(aes(y=V5,color="SDE+sd")) +geom_line(aes(y=V6,color="SDE-sd"))+geom_line(aes(y=V7,color="MSM"))+geom_line(aes(y=V8,color="MSM+sd"))+geom_line(aes(y=V9,color="MSM-sd")) + scale_shape_discrete(solid=T)   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x="Time",y=expression("Number of individuals in "*"S"[2])) 

ggsave(paste(currpath,"/n2.tiff",sep="") , units="in", width=6, height=5, dpi=1200, compression = 'lzw')

plot.n3<-ggplot(data.n3, aes(V1,y=value,color=variable)) +  geom_line(aes(y=V2,color="ODE")) + geom_line(aes(y=V3,color="CM"))+geom_line(aes(y=V4,color="SDE"))+geom_line(aes(y=V5,color="SDE+sd")) +geom_line(aes(y=V6,color="SDE-sd"))+geom_line(aes(y=V7,color="MSM")) +geom_line(aes(y=V8,color="MSM+sd"))+geom_line(aes(y=V9,color="MSM-sd"))+ scale_shape_discrete(solid=T)   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x="Time",y=expression("Number of individuals in "*"S"[3])) 

ggsave(paste(currpath,"/n3.tiff",sep="") , units="in", width=6, height=5, dpi=1200, compression = 'lzw')

plot.n4<-ggplot(data.n4, aes(V1,y=value,color=variable)) +  geom_line(aes(y=V2,color="ODE")) + geom_line(aes(y=V3,color="CM"))+geom_line(aes(y=V4,color="SDE"))+geom_line(aes(y=V5,color="SDE+sd")) +geom_line(aes(y=V6,color="SDE-sd"))+geom_line(aes(y=V7,color="MSM")) +geom_line(aes(y=V8,color="MSM+sd"))+geom_line(aes(y=V9,color="MSM-sd"))+ scale_shape_discrete(solid=T)   + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank(), axis.line = element_line(colour = "black"))+labs(x="Time",y=expression("Number of individuals in "*"S"[4])) 

ggsave(paste(currpath,"/n4.tiff",sep="") , units="in", width=6, height=5, dpi=1200, compression = 'lzw')

multiplot(plot.n1, plot.n2, plot.n3,plot.n4,  cols=2)

##############################################################################
##############################################################################
#end of code
##############################################################################