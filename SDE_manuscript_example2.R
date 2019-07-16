##############################################################################
##############################################################################
#This code is to conduct analyses for the example in the manuscript titled
#"Adding noise to Markov cohort models"
#under review at MDM journal
#author: Rowan Iskandar
#email: rowan_iskandar@brown.edu
#last edited: July 8 2019
#The code for the numerical exercise is available under a GNU GPL license and can be found at https://github.com/rowaniskandar/SDE.
##############################################################################
##############################################################################
#loading relevant libraries
library(yuima) #critical for stochastic differential equation method
library(deSolve)
library(ggplot2)
library(reshape2)
library(expm)
library(pracma)
library(markovchain) #critical for microsimulation
library(matrixcalc)
library(openxlsx)
#setting directory name
currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
set.seed(12345)
##############################################################################
##############################################################################
##############################################################################
#parameter values for breast cancer toxicity model
#reference: Alarid‐Escudero, Fernando, Anne H. Blaes, and Karen M. Kuntz. "Trade‐offs Between Efficacy and Cardiac Toxicity of Adjuvant Chemotherapy in Early‐Stage Breast Cancer Patients: Do Competing Risks Matter?." The breast journal 23.4 (2017): 401-409.
##############################################################################
#basic dynamic parameters (time duration, population size)
popsize <- 1000 #population size
ns <- 7 #number of states
nMC <-10000 # number of micro sim iterations
t0 <- 0 #initial time
tfin <-50 #final time
dt <- 1 #time step
tpoints <- seq(t0,tfin,dt) #number of discrete time points
muDieASR = c(0.004808,0.005185,0.005562,0.005936,0.006323,0.006739,0.007208,0.007742,0.008349,0.009024,0.009754,0.010550,0.011452,0.012489,0.013714)
muDieASR = 1/30
muMets = 0.06
muTox = 0
muTox_ch = 0.0054769
muTox_ac = 0.0054769*2.8659
muDieMets = 0.282938767
muDieTox = 0.156064775
muDieASR = 1/30
rrC_noch = 1
rrC_ch = 0.216
rrC_ac = 0.23
#transition rates
c12=muTox_ac
c12_noch=0
c13=muMets
c13_ac=muMets*rrC_ac
c17=muDieASR
c24=muMets
c24_ac=muMets*rrC_ac
c25=muDieTox
c27=muDieASR
c34=muTox_ac
c34_noch=0
c36=muDieMets
c37=muDieASR
c47=muDieASR
c45=muDieTox
c46=muDieMets
#setting up the c vector of transition rates
params_noch=c(c12_noch,c13,c17,c24,c25,c27,c34_noch,c36,c37,c45,c46,c47)
params_ac=c(c12,c13_ac,c17,c24_ac,c25,c27,c34,c36,c37,c45,c46,c47)
#setting up q-matrix (used in microsimulation)
Q_noch=matrix(c(-(c12_noch+c13+c17), c12_noch, c13, 0, 0, 0, c17,
        0, -(c24+c25+c27), 0, c24, c25, 0, c27,
        0, 0, -(c34_noch+c36+c37), c34_noch, 0, c36, c37,
        0, 0, 0, -(c45+c46+c47), c45, c46, c47,
        0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0),nrow=ns,ncol=ns,byrow=TRUE)
Q_ac=matrix(c(-(c12+c13_ac+c17), c12, c13_ac, 0, 0, 0, c17,
              0, -(c24_ac+c25+c27), 0, c24_ac, c25, 0, c27,
              0, 0, -(c34+c36+c37), c34, 0, c36, c37,
              0, 0, 0, -(c45+c46+c47), c45, c46, c47,
              0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0),nrow=ns,ncol=ns,byrow=TRUE)
init <-c(popsize,0,0,0,0,0,0) #initial state configuration
##############################################################################
##############################################################################
##stochastic differential equation method
##use yuima package method
y<-1 #step size for Euler Maruyama
tpoints <- seq(t0,tfin,dt/y) #number of discrete time points
sol <- c("n1","n2","n3","n4","n5","n6","n7") #health states
#for correlated Wiener (not needed):
# rho <-0.9
# sigma <- matrix(c(1,rho,rho,rho,rho,rho,rho,
#                   rho,1,rho,rho,rho,rho,rho,
#                   rho,rho,1,rho,rho,rho,rho,
#                   rho,rho,rho,1,rho,rho,rho,
#                   rho,rho,rho,rho,1,rho,rho,
#                   rho,rho,rho,rho,rho,1,rho,
#                   rho,rho,rho,rho,rho,rho,1,
# ),nrow=ns,ncol=n2,byrow=TRUE)
# C <- chol(sigma)
#ac/chemo arm
#drift vector
A_ac <-c("-n1*(c12+c13_ac+c17)",
      "n1*c12-n2*(c24_ac+c25+c27)",
      "n1*c13_ac-n3*(c34+c36+c37)",
      "n2*c24_ac+n3*c34-n4*(c45+c46+c47)",
      "n2*c25+n4*c45",
      "n3*c36+n4*c46",
      "n1*c17+n2*c27+n3*c37+n4*c47")
#diffusion matrix
B_ac <-matrix(c("n1*(c12+c13_ac+c17)","-n1*c12","-n1*c13_ac",0 ,0,0 ,"-n1*c17",
           "-n1*c12","n1*c12+n2*(c24_ac+c25+c27)",0, "-n2*c24_ac", "-n2*c25",0 ,"-n2*c27",
           "-n1*c13_ac",0,"n1*c13+n3*(c34+c36+c37)", "-n3*c34", 0, "-n3*c36", "-n3*c37",
              0,"-n2*c24_ac", "-n3*c34", "n2*c24_ac+n3*c34+n4*(c45+c46+c47)", "-n4*c45","-n4*c46","-n4*c47",
              0,"-n2*c25", 0 , "-n4*c45", "n2*c25+n4*c45", 0, 0 ,
              0, 0, "-n3*c36","-n4*c46", 0 , "n3*c36+n4*c46" , 0,
           "-n1*c17","-n2*c27","-n3*c37","-n4*c47",0, 0, "n1*c17+n2*c27+n3*c37+n4*c47")
              ,nrow=7,ncol=7)
#no chemo arm
#drift vector
A_noch <-c("-n1*(c12_noch+c13+c17)",
         "n1*c12_noch-n2*(c24+c25+c27)",
         "n1*c13-n3*(c34+c36+c37)",
         "n2*c24+n3*c34-n4*(c45+c46+c47)",
         "n2*c25+n4*c45",
         "n3*c36+n4*c46",
         "n1*c17+n2*c27+n3*c37+n4*c47")
#diffusion matrix
B_noch <-matrix(c("n1*(c12_noch+c13+c17)","-n1*c12_noch","-n1*c13",0 ,0,0 ,"-n1*c17",
                "-n1*c12_noch","n1*c12_noch+n2*(c24+c25+c27)",0, "-n2*c24", "-n2*c25",0 ,"-n2*c27",
                "-n1*c13",0,"n1*c13+n3*(c34+c36+c37)", "-n3*c34", 0, "-n3*c36", "-n3*c37",
                0,"-n2*c24", "-n3*c34", "n2*c24+n3*c34+n4*(c45+c46+c47)", "-n4*c45","-n4*c46","-n4*c47",
                0,"-n2*c25", 0 , "-n4*c45", "n2*c25+n4*c45", 0, 0 ,
                0, 0, "-n3*c36","-n4*c46", 0 , "n3*c36+n4*c46" , 0,
                "-n1*c17","-n2*c27","-n3*c37","-n4*c47",0, 0, "n1*c17+n2*c27+n3*c37+n4*c47")
              ,nrow=7,ncol=7)

MCtrace_all <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
MCtrace_all2 <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
tic()
#sample path loop
for(m in 1:nMC){
  mod <- setModel(drift=A_ac,diffusion=B_ac,state.variable=sol,solve.variable = sol)
  samp <- setSampling(Initial=t0,Terminal=tfin, n=length(tpoints)-1)
  yui <-setYuima(model=mod,sampling=samp)
  sde.out <- simulate(yui,xinit=init,true.parameter=params_ac)
  sde<-get.zoo.data(sde.out)
  MCtrace_all[m, ,1]<-as.vector(sde$`Series 1`)
  MCtrace_all[m, ,2]<-as.vector(sde$`Series 2`)
  MCtrace_all[m, ,3]<-as.vector(sde$`Series 3`)
  MCtrace_all[m, ,4]<-as.vector(sde$`Series 4`)
  MCtrace_all[m, ,5]<-as.vector(sde$`Series 5`)
  MCtrace_all[m, ,6]<-as.vector(sde$`Series 6`)
  MCtrace_all[m, ,7]<-as.vector(sde$`Series 7`)
}
toc()
#data.sde<-as.data.frame(sde)
data.MCtrace_all <- as.data.frame(MCtrace_all)
dum1 <-colSums(MCtrace_all)
#calculate life expectancy
LEMC <- rep(0,nMC)
LEsum <- 0
for (i in 1:nMC){
  LEsum <- 0
  LEsum <- sum(MCtrace_all[i, ,1])+sum(MCtrace_all[i, ,2])+sum(MCtrace_all[i, ,3])+sum(MCtrace_all[i, ,4]) #only considers counts in non-death states
  LEMC[i] <- LEsum/(popsize*y)
}  
LEMC_mean <- mean(LEMC) #mean of life expectancy
LEMC_std <-std(LEMC) #standard deviation of life expectancy
####
meanMC <-dum1/nMC
varMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdposMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)
sdnegMC <- matrix(rep(0,length(tpoints)*ns), nrow=length(tpoints), ncol=ns)

for (m in 1:nMC){
  MCtrace_all2[m, ,1] <-  (MCtrace_all[m, ,1]-meanMC[,1])^2
  MCtrace_all2[m, ,2] <-  (MCtrace_all[m, ,2]-meanMC[,2])^2
  MCtrace_all2[m, ,3] <-  (MCtrace_all[m, ,3]-meanMC[,3])^2
  MCtrace_all2[m, ,4] <-  (MCtrace_all[m, ,4]-meanMC[,4])^2
  MCtrace_all2[m, ,5] <-  (MCtrace_all[m, ,5]-meanMC[,5])^2
  MCtrace_all2[m, ,6] <-  (MCtrace_all[m, ,6]-meanMC[,6])^2
  MCtrace_all2[m, ,7] <-  (MCtrace_all[m, ,7]-meanMC[,7])^2
}

dum2 <-colSums(MCtrace_all2)
varMC <-dum2/nMC
sdposMC <- meanMC

sdposMC[,1] <- meanMC[,1] + sqrt(varMC[,1])/sqrt(nMC)
sdnegMC[,1] <- meanMC[,1] - sqrt(varMC[,1])/sqrt(nMC)
sdposMC[,2] <- meanMC[,2] + sqrt(varMC[,2])/sqrt(nMC)
sdnegMC[,2] <- meanMC[,2] - sqrt(varMC[,2])/sqrt(nMC)
sdposMC[,3] <- meanMC[,3] + sqrt(varMC[,3])/sqrt(nMC)
sdnegMC[,3] <- meanMC[,3] - sqrt(varMC[,3])/sqrt(nMC)
sdposMC[,4] <- meanMC[,4] + sqrt(varMC[,4])/sqrt(nMC)
sdnegMC[,4] <- meanMC[,4] - sqrt(varMC[,4])/sqrt(nMC)
sdposMC[,5] <- meanMC[,5] + sqrt(varMC[,5])/sqrt(nMC)
sdnegMC[,5] <- meanMC[,5] - sqrt(varMC[,5])/sqrt(nMC)
sdposMC[,6] <- meanMC[,6] + sqrt(varMC[,6])/sqrt(nMC)
sdnegMC[,6] <- meanMC[,6] - sqrt(varMC[,6])/sqrt(nMC)
sdposMC[,7] <- meanMC[,7] + sqrt(varMC[,7])/sqrt(nMC)
sdnegMC[,7] <- meanMC[,7] - sqrt(varMC[,7])/sqrt(nMC)

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
index_pull<-seq(1,tfin*y+y,dt*y)
data.SDEyui.mean <- data.frame(meanMC[index_pull,])
data.SDEyui.var <- data.frame(varMC[index_pull,])
data.SDEyui.sdpos <- data.frame(sdposMC[index_pull,])
data.SDEyui.sdneg <- data.frame(sdnegMC[index_pull,])
##############################################################################
##############################################################################
#microsimulation approach
##############################################################################
tpoints <- seq(t0,tfin,dt)
sickStates <- c("H","T","M","TM","DT","DTM","DO") #health state labels
Q=Q_ac
#transition probability matrix, calculated from the Q matrix
trprob <- expm(Q*tstep)
P <- trprob
#check stochastic matrix, make sure row sum = 1
last_index = rep(0,ns)
row_sum = rep(0,ns)
cum_sum = matrix(rep(0,ns*ns),nrow=ns,ncol=ns)
for (i in 1:ns){
  row_sum[i]=0
  for (j in 1:ns){
    if (P[i,j]!=0 & row_sum[i] < 1){
    row_sum[i]=row_sum[i]+P[i,j]
    cum_sum[i,j]=row_sum[i]
    last_index[i]=j
    }
    else{
      row_sum[i]=row_sum[i]
    }
  }
}
for (i in 1:ns){
  if (sum(P[i,]) !=1){
  P[i,last_index[i]]=1-cum_sum[i,last_index[i]-1]
  }
}
trprob <- P
sickP <-matrix(data = trprob,  nrow = ns, ncol=ns, dimnames = list(sickStates,sickStates))
#use markovchain package
sickMC <- new("markovchain", states = sickStates , byrow=TRUE, transitionMatrix =sickP , name = "sick")

N <- popsize
tfinal <- length(tpoints)
initialState <- c(1,0,0,0,0,0,0)

SHtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SHbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
STtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
STbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SMtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SMbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
STMtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
STMbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDTtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDTbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDMtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDMbool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDOtrace <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
SDObool <- matrix(rep(0,N*length(tpoints)), nrow = N, ncol = tfinal)
MCtrace <- matrix(rep(0,length(tpoints)*ns), nrow = length(tpoints), ncol = ns)
MCtrace_all <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))
MCtrace_all2 <-array(rep(0, nMC*length(tpoints)*ns), dim=c(nMC, length(tpoints), ns))

tic() #timer start
#start outer outer loop (sample path)
for (m in 1:nMC){
  #start individual path
  for (i in 1:N){
    sicktrace[i,] <- rmarkovchain(n = tfinal, object = sickMC, t0 = "H") #generate markov chain using markovchain package
  }
  for (i in 1:N){
    for (j in 1:tfinal){
      if (sicktrace[i,j]=="H"){SHbool[i,j]=1}
      else {SHbool[i,j]=0}
      if (sicktrace[i,j]=="T"){STbool[i,j]=1}
      else {STbool[i,j]=0}
      if (sicktrace[i,j]=="M"){SMbool[i,j]=1}
      else {SMbool[i,j]=0}
      if (sicktrace[i,j]=="TM"){STMbool[i,j]=1}
      else {STMbool[i,j]=0}
      if (sicktrace[i,j]=="DT"){SDTbool[i,j]=1}
      else {SDTbool[i,j]=0}
      if (sicktrace[i,j]=="DM"){SDMbool[i,j]=1}
      else {SDMbool[i,j]=0}
      if (sicktrace[i,j]=="DO"){SDObool[i,j]=1}
      else {SDObool[i,j]=0}
    }
  }
  MCtrace_all[m, ,1] <- colSums(SHbool, dims=1)
  MCtrace_all[m, ,2] <- colSums(STbool, dims=1)
  MCtrace_all[m, ,3] <- colSums(SMbool, dims=1)
  MCtrace_all[m, ,4] <- colSums(STMbool, dims=1)
  MCtrace_all[m, ,5] <- colSums(SDTbool, dims=1)
  MCtrace_all[m, ,6] <- colSums(SDMbool, dims=1)
  MCtrace_all[m, ,7] <- colSums(SDObool, dims=1)
}
toc()
dum1 <-colSums(MCtrace_all)
meanMC <-dum1/nMC

#calculate life expectancy
LEMC <- rep(0,nMC)
LEsum <- 0
for (i in 1:nMC){
  LEsum <- 0
  LEsum <- sum(MCtrace_all[i, ,1])+sum(MCtrace_all[i, ,2])+sum(MCtrace_all[i, ,3])+sum(MCtrace_all[i, ,4])
  LEMC[i] <- LEsum/popsize
}  
LEMSM_mean <- mean(LEMC)
LEMSM_std <- std(LEMC)

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
#output dataframes of relevant statistics
data.MC.mean <- data.frame(meanMC_d)
data.MC.var <- data.frame(varMC)
data.MC.sdpos <- data.frame(sdposMC_d)
data.MC.sdneg <- data.frame(sdnegMC_d)

##############################################################################
##############################################################################
#end of code
##############################################################################
