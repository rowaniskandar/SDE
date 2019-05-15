######################################################
#Example 1 in manuscript titled:
#Adding noise to Markov cohort models
#requires the following package/module
include("./CohMod.jl")
using Plots
# using GraphPlot
# using LightGraphs
#testing CohMod.jl module
#Rowan Iskandar (rowan.iskandar@gmail.com)
#start 10042019
#two examples
######################################################
#example 1: 4-state example in
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205543
#set up parameters
c12 = 0.05
c13 = 0.01
c14 = 0.001
c23 = 0.1
c24 = 0.05
c34 = 1
popsize = 1000
nMC = 10000 # number of micro sim iterations
t0 = 0 #initial time
tfin = 50 #final time
dt = 1 #time step
dt_scale = 1
c = [c12 ;c13 ;c14 ;c23 ;c24 ;c34 ]
pop_init = [popsize 0 0 0] #initial state configuration
#generator matrix
Q = [-(c12+c13+c14) c12 c14 c14;0 -(c23+c24) c23 c24;0 0 -c34 c34 ;0 0 0 0]
d=[-1 1 0 0;-1 0 1 0 ;-1 0 0 1 ;0 -1 1 0 ;0 -1 0 1;0 0 -1 1]
state_names=["S1" "S2" "S3" "S4"]
#call module and set global parameter values
CohMod.SetParams(state_names,t0,tfin,dt,dt_scale,Q,c,d,pop_init,popsize)
LEdims=[1 2 3]
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
#################################################################
#end of file
