#
# Persephone: Epidemics Spread in Hierarchical Network Models
# Copyright (C) 2020  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

using LightGraphs, Pandas, Plots, GraphPlot, Compose, Cairo, Fontconfig, Measures, Printf, Colors

using Adjacently.util
using Adjacently.io
using Adjacently.pr
using Adjacently.graph

include("./utils.jl")

##########
# simulation variables
##########

# whole period Copenhagen graph
NV = 692
NE = 79530

# 
# Reproduction number	 R0 	2.4
# Incubation period (days)	 τincubation 	5.1
# Infectious period (days)	 τinfectious 	3.3

# run 1: (20, 100, 0.8, 1:10 infected)
# run 2: (1469, 950, 0.4, 1:10 infected)
# run 3: (1469, 950, 0.4, highest PR infected)
# incubation period (5 min bins)
E_t = 1469
# infectious period (5 min bins)
I_t = 950
# probability of getting infected in 1 interaction (5 min bin)
P_i = 0.4

EPS = 1e-8

# S,E,I,R
VCOLOR = [colorant"blue" colorant"orange" colorant"red" colorant"grey"]

# save anim
SAVE_ANIM = false

N_BINS = 8000

##########
# initializations
##########

@info("loading bt-detected-internal-notime")
g_all = SimpleGraph{UInt32}()
oni = load_adjacency_list_from_csv(UInt32, g_all, "./data/bt-detected-internal-notime.csv")
# write_mgs3_graph(g_all, "bt-detected-internal-notime")

# initialize stat
stat = Stat(NV,0,0,0)

# initalize vstats
vstats = Dict{UInt32,VState}()
for v in 1:NV
	vstats[convert(UInt32,v)] = VState(S::Status,0)
end

# compute pagerank
P = get_sparse_P_matrix(g_all)
pr = zeros(Float64,NV)
pr = PR(P, epsilon=EPS)

# id: 251
val_pr_max, ind_pr_max = findmax(pr)
@info("PR max: ", val_pr_max, ind_pr_max)

# infect highest PR node
vstats[convert(UInt32,ind_pr_max)].status = I::Status
s2e(stat)
e2i(stat)

# infect 10 nodes
#for j in 1:10
#	vstats[convert(UInt32,j)].status = I::Status
#	s2e(stat)
#	e2i(stat)
#end

# set spring layout
l_x, l_y = spring_layout(g_all; C=80)

g_base = SimpleGraph{UInt32}()
add_vertices!(g_base,NV)

# time varies fom 0s to 2418900s to 300s increments (5 min) -> 8063 periods
# 692 vs / 79530 edges
# CSV file headers: source, target, time
@info("loading bt-detected-internal")
es = Pandas.read_csv("./data/bt-detected-internal.csv")

#if SAVE_ANIM
#    anim = Animation()
#end

## load edge time bins
#ebins = Dict{Int,Array{UInt32,2}}()
#
#@info("loading time bin edges")
#for i in 0:(N_BINS - 1)
#    ebins[i] = get_bin_edges(es, i)
#end

SA = zeros(Int, N_BINS)
EA = zeros(Int, N_BINS)
IA = zeros(Int, N_BINS)
RA = zeros(Int, N_BINS)

##########
# simulation
##########

for i in 0:(N_BINS - 1)
    # create temporal graph of current interactions
    bin_edges = get_bin_edges(UInt32,es,i)
    g = get_temporal_graph(UInt32,NV,oni,bin_edges) 

    #if SAVE_ANIM
    #	# enum are 0-based!?
    #	# membership = [Int(vstats[v].status)+1 for v in 1:NV]
    #	membership = [Int(vstats[v].status) for v in 1:NV]
    #	# membership color
    #	nodefillc = VCOLOR[membership]

    #	p = gplot(g_base, l_x, l_y, nodefillc=nodefillc)
    #	output = compose(p)
    #
    #	j = length(anim.frames) + 1
    #	tmpfilename = joinpath(anim.dir, @sprintf("%06d.png", j))
    #	Compose.draw(PNG(tmpfilename), output)
    #	push!(anim.frames, tmpfilename)
    #end
    
    SA[i+1] = stat.ns
    EA[i+1] = stat.ne
    IA[i+1] = stat.ni
    RA[i+1] = stat.nr

    # spread infection
    spread(g,i,vstats,stat,e_t=E_t,p_i=P_i,i_t=I_t)

    # garbage collect
    if i % 100 == 0
    	GC.gc()
    end
end

#if SAVE_ANIM
#	gif(anim, "./plots/animation.gif", fps = 10)
#end

Plots.plot(title = "SEIR Simulation", xlabel = "Time", ylabel = "# individuals")
Plots.plot([SA EA IA RA], label = ["S" "E" "I" "R"], linecolor = VCOLOR)

png("./plots/time-series.png")

