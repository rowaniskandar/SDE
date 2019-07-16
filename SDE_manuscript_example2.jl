######################################################
#Examples in manuscript titled:
#Adding noise to Markov cohort models
#include functions from the following package
include("./CohModjl")
using Plots
using CSV
using JLD2
#Rowan Iskandar (rowan.iskandar@gmail.com)
#start 10042019
#last revised: July 8 2019
######################################################
######################################################
######################################################
#example 2: bc toxicity model
#https://onlinelibrary.wiley.com/doi/abs/10.1111/tbj.12757
#STEP 1
#list of states (s=7)
#1: ned, 2: mets, 3: tox, 4: metstox, 5: dieoth, 6: dietox, 7: diemets
#enumerate transitions (r=12)
#from healthy
#1-2
#1-3
#1-7
#from tox
#2-4
#2-5
#2-7
#from mets
#3-4
#3-6
#3-7
#from mets+tox
#4-5
#4-6
#4-7
#create matrix of change (matrix d)
#row: transition type
#column: change to individual state: increase/decrease by one (+1,-1), or no change (0)
d=[-1 1 0 0 0 0 0;
-1 0 1 0 0 0 0;
-1 0 0 0 0 0 1;
0 -1 0 1 0 0 0;
0 -1 0 0 1 0 0;
0 -1 0 0 0 0 1;
0 0 -1 1 0 0 0;
0 0 -1 0 0 1 0;
0 0 -1 0 0 0 1;
0 0 0 -1 1 0 0;
0 0 0 -1 0 1 0;
0 0 0 -1 0 0 1]
#transition rates common to all intervention arms
muMets = 0.06
muTox = 0
muTox_ch = 0.0054769
muTox_ac = 0.0054769*2.8659
muDieMets = 0.282938767
muDieTox = 0.156064775
#getting lifetable
path = dirname(dirname(@__FILE__))
file_lifetable = string(path,"/codes/lifetable.csv")
lifetable_df= CSV.read(file_lifetable,datarow=2)
lifetable =convert(Array,lifetable_df)
#relative risks
rrC_noch = 1
rrC_ch = 0.216
rrC_ac = 0.23
#noch arm-specific transition rates
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
c45=muDieTox
c46=muDieMets
c47=muDieASR
#vector of transition rates (used in propensitive vector v in CohMod.jl)
c_noch=[c12_noch;c13;c17;c24;c25;c27;c34_noch;c36;c37;c45;c46;c47]
c_ac=[c12;c13_ac;c17;c24_ac;c25;c27;c34;c36;c37;c45;c46;c47]
#Q matrix or generator matrix (used in microsimulation)
Q_noch=[-(c12_noch+c13+c17) c12_noch c13 0 0 0 c17;
0 -(c24+c25+c27) 0 c24 c25 0 c27;
0 0 -(c34_noch+c36+c37) c34_noch 0 c36 c37;
0 0 0 -(c45+c46+c47) c45 c46 c47;
0 0 0 0 0 0 0;
0 0 0 0 0 0 0;
0 0 0 0 0 0 0]
Q_ac=[-(c12+c13_ac+c17) c12 c13_ac 0 0 0 c17;
0 -(c24_ac+c25+c27) 0 c24_ac c25 0 c27;
0 0 -(c34+c36+c37) c34 0 c36 c37;
0 0 0 -(c45+c46+c47) c45 c46 c47;
0 0 0 0 0 0 0;
0 0 0 0 0 0 0;
0 0 0 0 0 0 0]
#basic dynamic parameters (time duration, population size)
popsize = 1000  #population size
nMC = 10000 # number of monte carlo samples both for microsimulation and SDE
t0 = 0 #initial time
tfin = 45 #final time
dt = 1 #time step
dt_scale = 1 #for finer time steps in Euler Maruyama algorithm
pop_init = [popsize 0 0 0 0 0 0]  #initial state configuration
state_names=["NED" "Tox" "Mets" "Mets+Tox" "DieTox" "DieMets" "DieASR"]
LEdims=[1 2 3 4] #flag which column is non-death state (to calculate life expectancy or survival)
#call module and set global parameter values
#simulating ac arm
CohMod_ver2.SetParams(state_names,t0,tfin,dt,dt_scale,c12,c13_ac,c17,c24_ac,c25,c27,c34,c36,c37,c45,c46,c47,d,pop_init,popsize,lifetable)
#######################################################################
#call different methods
#first method: SDE
#####STEP 2-5 are contained in the file CohMod.jl
@time begin #calculating elapsed time
elapsedtime_SDE,LEmean_SDE, LEstd_SDE, tpoints, OSmean_SDE, OSstd_SDE, SDEmean_noch, SDEnegsd_noch, SDEpossd_noch = CohMod_ver2.StochDiffEquation(nMC,LEdims)
#tpoints,  SDEmean_noch, SDEnegsd_noch, SDEpossd_noch = CohMod.StochDiffEquation(nMC,LE,LEdims)
end
#######################################################################
#second method: MicroSimulation
@time begin
elapsedtime_MC, LEmean_MSM, LEstd_MSM, tpoints, OSmean_MSM, OSstd_MSM, MCmean_noch, MCnegsd_noch, MCpossd_noch   = CohMod_ver2.MicroSimulation(nMC,LEdims)
end
results_1000_10000 = [LEmean_MSM  LEstd_MSM elapsedtime_MC; LEmean_SDE LEstd_SDE elapsedtime_SDE; 0 popsize nMC]
@save "tmpfile_1000_10000_2.jld" results_1000_10000
#################################################################
#end of file
