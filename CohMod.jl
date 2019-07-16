module CohMod
###############################################################
# Copyright 2019 Rowan Iskandar
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#################################################################
#CohMod module include the following modeling methods to represent
# the trajectories of a cohort of patients (Markov state-transition cohort model)
#used in
#1. differential equation (ODE) - using julia package DifferentialEquations
#2. stochastic differential equation (SDE)
#3. microsimulation
#4. master equation - gillespie
#5. poisson representation
#The relationship among the methods are discussed in
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205543
#and the forthcoming
#"Adding noise to Markov state-transition cohort models"
##################################################################
#use packages
using DifferentialEquations, DiffEqBiological
using LinearAlgebra, Statistics, Compat
using Distributions #, Printf,Random
#using LightGraphs, MetaGraphs, GraphPlot, Colors, Compose
# using Fontconfig, Cairo
using PoissonRandom
using Dates
using Random
using Nemo
#SetParams initialize module with required global parameters
#t0:pop_initial time, tfin:final time, dt:time step
#dt_scale:scale dt for SDE
#c: vector of rate constants, dim(c): number of transitions
#pop_init: vector of pop_initial state-config, dim(pop_init): number of states
#nMC: number of monte-carlo simulations
#d: matrix of stoichiometries, dim(d): dim(c)x s
#N: N(t) current state-configuration, or N(t)
#Q: generator (Q) matrix which contains elements of c, must be square
#################################################################
#this function set the basic parameters + method-specific parameters
function SetParams(state_names_in,t0_in,tfin_in,dt_in,dt_scale_in,c12_in,c13_ac_in,c17_in,c24_ac_in,c25_in,c27_in,c34_in,c36_in,c37_in,c45_in,c46_in,c47_in,d_in,pop_init_in,popsize_in,lifetable_in)
     global state_names=state_names_in
     global t0 = t0_in
     global tfin=tfin_in
     global dt=dt_in
     global dt_scale=dt_scale_in
     #global c=c_in
     global d=d_in
     global pop_init=pop_init_in
     global popsize=popsize_in
     #global nMC=nMC_in
     #global Q=Q_in
     global c12=c12_in
     #c12_noch=0
     #c13=c13_in
     global c13_ac=c13_ac_in
     global c17=c17_in
     #c24=c24_in
     global c24_ac=c24_ac_in
     global c25=c25_in
     global c27=c27_in
     global c34=c34_in
     #global c34_noch=c34_noch_in
     global c36=c36_in
     global c37=c37_in
     global c45=c45_in
     global c46=c46_in
     global c47=c47_in
     global num_s = length(pop_init)
     global num_r = size(d,1)
     global tpoints = collect(t0:dt:tfin)
     global tpoints_SDE= collect(t0:dt/dt_scale:tfin)
     global num_t = length(tpoints)
     global num_t_SDE = length(tpoints_SDE)
     global index_pull=collect(1:dt*dt_scale:tfin*dt_scale+dt_scale-1)
     global lifetable=lifetable_in
end
#################################################################
#stochastic differential equation (SDE)
#arbirary num_s (number of states), num_r (number of transitions),
#arbitrary transition rates, d (matrix of changes), P (transition probability matrix)
#start 11042019
#calculate propensities of reactions at time t (depends on N)
#in the manuscrip this is the vector v which a function of c, d, and N)
function vprop(N,c,d)
    v=zeros(num_r)
    for i=1:num_r
        index_n=1 #for mono-reaction
        index_find=false
        while index_find==false
            if d[i,index_n]<0
                index_find=true
                v[i]=c[i]*N[index_n]
            else
                index_n=index_n+1
            end
        end
    end
    return v
end

#calculate mean change (drift vector)
function change_mean(N,c,d)
    A=vprop(N,c,d)'*d
    return A
end
#calculate variance of change (diffuion matrix)
function change_variance(N,c,d)
    B=zeros(Float16,num_s,num_s)
    for i=1:num_s
        for j=1:num_s
            sum=0
            for k=1:num_r
                sum=sum+vprop(N,c,d)[k]*d[k,i]*d[k,j]
            end
            B[i,j]=sum
        end
    end
    return B
end
#solving SDE using Euler Maruyama method
function StochDiffEquation(nMC,LEdims)
    Random.seed!(1234)
    start=time()
    pop_trace = zeros(Float64,num_s,num_t_SDE)
    MC_pop_trace = zeros(Float64,nMC, num_s,num_t_SDE)
    MC_pop_mean = zeros(Float64,num_s,num_t_SDE)
    MC_pop_std = zeros(Float64,num_s,num_t_SDE)
    MC_LE_mean = zeros(Float64,num_s,num_t_SDE)
    MC_LE_std = zeros(Float64,num_s,num_t_SDE)
    OS=zeros(Float64,nMC,num_t)
    OSmean=zeros(Float64,num_t)
    OSstd=zeros(Float64,num_t)
    LEsum=zeros(Float64,nMC)
    LEmean::Float64=0
    LEstd::Float64=0
    pop_trace[:,1] = pop_init
    dist_norm = Normal(0,1)
    #monte carlo samples (sample paths)
    for n=1:nMC
        for t=2:num_t_SDE
            dW=rand(dist_norm,num_s)
            #update time dependent parameters
            c17=c27=c37=c47=lifetable[t-1,2]
            c=[c12;c13_ac;c17;c24_ac;c25;c27;c34;c36;c37;c45;c46;c47]
            #calculate drift
            A=change_mean(pop_trace[:,t-1],c,d)
            #calculate diff
            B=change_variance(pop_trace[:,t-1],c,d)
            Beigvals=eigvals(B)
            for m=1:num_s
                if Beigvals[m]<0
                    Beigvals[m]=0
                end
            end
            Bdiag=Diagonal(Beigvals)
            L=eigvecs(B)
            Linv=inv(L)
            Bsqrt=L*sqrt(Bdiag)*Linv
            pop_trace[:,t]=pop_trace[:,t-1]+A'*(dt/dt_scale)+Bsqrt*sqrt(dt/dt_scale)*dW
            #avoiding negative population
            #current approach (10/05/2019)
            for j=1:num_s
                if pop_trace[j,t]<0
                    pop_trace[j,t]=0
                end
            end
        end
        MC_pop_trace[n,:,:]=pop_trace

        for i=1:length(LEdims)
            LEsum[n]=LEsum[n]+sum(pop_trace[LEdims[i],:])
        end
        LEsum[n]=LEsum[n]/popsize
    end
    for i=1:num_s
        x=MC_pop_trace[:,i,:]
        MC_pop_mean[i,:]=mean(x;dims=1)
        MC_pop_std[i,:]=std(x;dims=1)
    end
    #return tpoints, MC_pop_mean, MC_pop_mean+MC_pop_std, MC_pop_mean-MC_pop_std
    #return tpoints, MC_pop_mean[:,index_pull], MC_pop_mean[:,index_pull]+MC_pop_std[:,index_pull], MC_pop_mean[:,index_pull]-MC_pop_std[:,index_pull]
    #calculate overall survival
    for j=1:nMC
        for i=1:length(LEdims)
            OS[j,:]=OS[j,:]+MC_pop_trace[j,LEdims[i],:]
        end
    end
    OSmean=mean(OS/popsize;dims=1)
    OSstd=std(OS/popsize;dims=1)
    #calculate life expectancy
    LEmean=mean(LEsum)
    LEstd=std(LEsum)
    #####
    elapsedtime=time()-start
    return elapsedtime, LEmean,  LEstd, tpoints, OSmean, OSstd, MC_pop_mean[:,index_pull], MC_pop_mean[:,index_pull]+MC_pop_std[:,index_pull], MC_pop_mean[:,index_pull]-MC_pop_std[:,index_pull]

end
#################################################################
#master equation - gillespie or SSA method
#arbirary num_s, num_r, c, d, P
#use DiffEqBiological.jl API
#output tpoints, mean, mean-std, mean+std
function MasterEqGil(nMC)
    #native
    pop_trace = zeros(Int32,1,num_s+1) #extra column to store time
    MC_pop_trace = zeros(Int32,nMC, num_s,num_t)
    MC_pop_mean = zeros(Float16,num_s,num_t)
    MC_pop_std = zeros(Float16,num_s,num_t)
    t=0
    t_next=0
    index_reaction=1
    pop_trace[1,:]=[t pop_init]
    pop_now = pop_init
    while t<tfin
        vprop_now=vprop(pop_now,c,d)
        vprop_sum=sum(vprop_now)
        # dist_exp=Exponential(vprop_sum)
        # t_next=rand(dist_exp,1)
        U=rand()
        t_next=-log(1-U)/vprop_sum
        Pnow=vprop_now/vprop_sum
        print(Pnow)
        dists=Categorical(Pnow)
        index_k=rand(dists,1)
        pop_now=pop_now+d[index_reaction,:]'
        t=t+t_next[1]
        pop_trace=vcat(pop_trace,[t pop_now])
    end
    return pop_trace
end

# #################################################################
#microsimulation
#arbirary num_s, num_r, c, d, P
#output mean, mean-std, mean+std
#start 11042019
function MicroSimulation(nMC, LEdims)
    start=time()
    # setup the simulation
    MC_pop_trace = zeros(Int32, nMC,popsize,num_t)
    pop_trace = zeros(Int32,popsize,num_t)
    MC_pop_count = zeros(Int32,nMC, num_s,num_t)
    MC_pop_mean = zeros(Float16,num_s,num_t)
    MC_pop_std = zeros(Float16,num_s,num_t)
    MC_LE_mean = zeros(Float64,num_s,num_t_SDE)
    MC_LE_std = zeros(Float64,num_s,num_t_SDE)
    LEsum=zeros(Float64,nMC)
    OS=zeros(Float64,nMC,num_t)
    OSmean=zeros(Float64,num_t)
    OSstd=zeros(Float64,num_t)
    LEmean::Float64=0
    LEstd::Float64=0
    pop_count = zeros(Int32, num_s,num_t)
    X = fill(0, num_t) # allocate memory, or zeros(Int64, sample_size)
    for n=1:nMC
        for i=1:popsize
            X[1] = 1 # set the pop_initial state
            for t in 2:num_t
                #update time-dependent parameters
                c17=c27=c37=c47=lifetable[t-1,2]
                Q=[-(c12+c13_ac+c17) c12 c13_ac 0 0 0 c17;
                0 -(c24_ac+c25+c27) 0 c24_ac c25 0 c27;
                0 0 -(c34+c36+c37) c34 0 c36 c37;
                0 0 0 -(c45+c46+c47) c45 c46 c47;
                0 0 0 0 0 0 0;
                0 0 0 0 0 0 0;
                0 0 0 0 0 0 0]
                P=exp(Q)
                #check stochastic matrix, make sure row sum = 1
                last_index = zeros(Int8,num_s)
                row_sum = zeros(num_s)
                cum_sum = zeros(num_s,num_s)
                for k=1:num_s
                    row_sum[k]=0
                    for j=1:num_s
                        if P[k,j]!==0 && row_sum[k] < 1
                            row_sum[k]=row_sum[k]+P[k,j]
                            cum_sum[k,j]=row_sum[k]
                            last_index[k]=j
                        else
                            row_sum[k]=row_sum[k]
                        end
                    end
                end
                for k=1:num_s
                    if sum(P[k,:]) !==1
                        P[k,last_index[k]]=1-cum_sum[k,last_index[k]-1]
                    end
                end
                @assert size(P)[1] == size(P)[2] # square required
                N = size(P)[1] # should be square
                # create vector of discrete RVs for each row
                dists = [Categorical(P[i, :]) for i in 1:N]
                dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
                X[t] = rand(dist) # draw new value
            end
            pop_trace[i,:] = X
        end
        # collecting results: counting number of people in each state for each stoichiometries
        MC_pop_trace[n,:,:] = pop_trace
        for t=1:tfin
            for i=1:num_s
                pop_count[i,t] = count(pop_trace[:,t].== i)
                #pop_count[i,t] = count(pop_trace[pop_trace[:,t].== i,t])
            end
        end
        MC_pop_count[n,:,:]=pop_count

        for i=1:length(LEdims)
            LEsum[n]=LEsum[n]+sum(pop_count[LEdims[i],:])
        end
        LEsum[n]=LEsum[n]/(popsize*dt_scale)
    end
    for i=1:num_s
        x=MC_pop_count[:,i,:]
        MC_pop_mean[i,:]=mean(x;dims=1)
        MC_pop_std[i,:]=std(x;dims=1)
    end
    #calculate overall survival
    for j=1:nMC
        for i=1:length(LEdims)
            OS[j,:]=OS[j,:]+MC_pop_count[j,LEdims[i],:]
        end
    end
    OSmean=mean(OS/popsize;dims=1)
    OSstd=std(OS/popsize;dims=1)
    #calculate life expectancy
    LEmean=mean(LEsum)
    LEstd=std(LEsum)
    #####

    elapsedtime = time()-start
    return elapsedtime, LEmean, LEstd, tpoints, OSmean, OSstd, MC_pop_mean, MC_pop_mean+MC_pop_std, MC_pop_mean-MC_pop_std

end

#################################################################
#end of module
end
