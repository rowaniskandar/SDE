######################################################
#Examples in manuscript titled:
#Adding noise to Markov cohort models
#include functions from the following package
include("./CohMod.jl")
using Plots
# using GraphPlot
# using LightGraphs
#testing CohMod.jl module
#Rowan Iskandar (rowan.iskandar@gmail.com)
#start 10042019
######################################################
######################################################
#example 2: bc toxicity model
#https://onlinelibrary.wiley.com/doi/abs/10.1111/tbj.12757
#list of states (num_s=7)
#1: ned, 2: mets, 3: tox, 4: metstox, 5: dieoth, 6: dietox, 7: diemets
#enumerate transitions (num_r=12)
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
#noch arm
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

c_noch=[c12_noch;c13;c17;c24;c25;c27;c34_noch;c36;c37;c45;c46;c47]
c_ac=[c12;c13_ac;c17;c24_ac;c25;c27;c34;c36;c37;c45;c46;c47]
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
popsize = 1000
nMC = 100000 # number of micro sim iterations
t0 = 0 #initial time
tfin = 45 #final time
dt = 1 #time step
dt_scale = 1
pop_init = [popsize 0 0 0 0 0 0] #initial state configuration
state_names=["NED" "Tox" "Mets" "Mets+Tox" "DieTox" "DieMets" "DieASR"]
LEdims=[1 2 3 4]
#call module and set global parameter values
#noch arm
CohMod.SetParams(state_names,t0,tfin,dt,dt_scale,Q_ac,c_ac,d,pop_init,popsize)
# g, adjacency_matrix, edgelabel=CohMod.SetGraph()
# gplot(g,nodelabel=state_names,edgelabel=edgelabel)
#call different methods
@time begin
elapsedtime_SDE,LEmean_SDE, LEstd_SDE, tpoints, OSmean_SDE, OSstd_SDE, SDEmean_noch, SDEnegsd_noch, SDEpossd_noch = CohMod.StochDiffEquation(nMC,LEdims)
#tpoints,  SDEmean_noch, SDEnegsd_noch, SDEpossd_noch = CohMod.StochDiffEquation(nMC,LE,LEdims)
end
# plot(tpoints,SDEmean_noch[1,:])
# plot!(tpoints,SDEmean_noch[2,:])
# plot!(tpoints,SDEmean_noch[3,:])
# plot!(tpoints,SDEmean_noch[4,:])
# plot!(tpoints,SDEmean_noch[5,:])
# plot!(tpoints,SDEmean_noch[6,:])
# plot!(tpoints,SDEmean_noch[7,:])
@time begin
Ptrans, elapsedtime_MC, LEmean_MSM, LEstd_MSM, tpoints, OSmean_MSM, OSstd_MSM, MCmean_noch, MCnegsd_noch, MCpossd_noch   = CohMod.MicroSimulation(nMC,LEdims)
#tpoints, MCmean_noch, MCnegsd_noch, MCpossd_noch  = CohMod.MicroSimulation(nMC)
end
plot(tpoints[1:50],OSmean_MSM[1:50])
plot!(tpoints[1:50],OSmean_SDE[1:50])
plot!(tpoints[1:50],OSmean_SDE[1:50]+OSstd_SDE[1:50])
plot!(tpoints[1:50],OSmean_SDE[1:50]-OSstd_SDE[1:50])

plot(tpoints[1:50],MCmean_noch[1,1:50])
plot!(tpoints[1:50],SDEmean_noch[1,1:50])
plot!(tpoints[1:50],MCmean_noch[2,1:50])
plot!(tpoints[1:50],MCmean_noch[3,1:50])
plot!(tpoints[1:50],MCmean_noch[4,1:50])
plot!(tpoints[1:50],MCmean_noch[5,1:50])
plot!(tpoints[1:50],MCmean_noch[6,1:50])
plot!(tpoints[1:50],MCmean_noch[7,1:50])
@time begin
ODE_noch, ODEstiff_noch = CohMod.DiffEquation(, , Tsit5())
end
plot(ODE_noch)
plot(ODEstiff_noch)
#ac arm
CohMod.SetParams(state_names,t0,tfin,dt,dt_scale,Q_ac,c_ac,d,pop_init,popsize)
#using DiffEqBiological
# rn = @reaction_network begin
#     c12, S1 --> S2
#     c13, S1 --> S3
#     c15, S1 --> S5
#     c24, S2 --> S4
#     c25, S2 --> S5
#     c27, S2 --> S7
#     c34, S3 --> S4
#     c35, S3 --> S5
#     c36, S3 --> S6
#     c45, S4 --> S5
#     c46, S4 --> S6
#     c47, S4 --> S7
# end
# jumprateexpr=rn.jump_rate_expr
# jumpexpr=rn.jump_affect_expr
# jumpode=rn.odefun
# trans=rxtospecies_depgraph(rn)
# specious=speciestorx_depgraph(rn)
# rx=rxtorx_depgraph(rn)
#call different methods
@time begin
tpoints,  SDEmean_ac, SDEnegsd_ac, SDEpossd_ac = CohMod.StochDiffEquation(nMC)
end
plot(tpoints,SDEmean_ac[1,:])
plot!(tpoints,SDEmean_ac[2,:])
plot!(tpoints,SDEmean_ac[3,:])
plot!(tpoints,SDEmean_ac[4,:])
plot!(tpoints,SDEmean_ac[5,:])
plot!(tpoints,SDEmean_ac[6,:])
plot!(tpoints,SDEmean_ac[7,:])

@time begin
tpoints, MCmean_ac, MCnegsd_ac, MCpossd_ac  = CohMod.MicroSimulation(nMC)
end
plot(tpoints[1:50],MCmean_ac[1,1:50])
plot!(tpoints[1:50],MCmean_ac[2,1:50])
plot!(tpoints[1:50],MCmean_ac[3,1:50])
plot!(tpoints[1:50],MCmean_ac[4,1:50])
plot!(tpoints[1:50],MCmean_ac[5,1:50])
plot!(tpoints[1:50],MCmean_ac[6,1:50])
plot!(tpoints[1:50],MCmean_ac[7,1:50])
@time begin
ODE_ac, ODEstiff_ac = CohMod.DiffEquation(, , Tsit5())
end
plot(ODE_ac)
plot(ODEstiff_ac)
#################################################################
#end of file
