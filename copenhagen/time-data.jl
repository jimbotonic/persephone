#
# Persephone: Epidemics Spread in Hierarchical Network Models
# Copyright (C) 2016-2020  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using LightGraphs, Pandas, Plots, GraphPlot, Compose, Cairo, Fontconfig, Measures, Printf, Colors, GR

using Adjacently.util
using Adjacently.io
using Adjacently.pr
using Adjacently.graph

# susceptible (S), exposed (E), infected (I), and resistant (R)
@enum Status S=1 E=2 I=3 R=4

# node state
mutable struct VState
	status::Status
	time::UInt32
end

# current stats
mutable struct Stat
	ns::Int
	ne::Int
	ni::Int
	nr::Int
end

function s2e(stat::Stat)
	if stat.ns > 0
		stat.ns = stat.ns - 1
		stat.ne = stat.ne + 1
	else
		error("# susceptible is 0!")
	end
end
function e2i(stat::Stat)
	if stat.ne > 0
		stat.ne = stat.ne - 1
		stat.ni = stat.ni + 1
	else
		error("# exposed is 0!")
	end
end
function i2r(stat::Stat)
	if stat.ni > 0
		stat.ni = stat.ni - 1
		stat.nr = stat.nr + 1
	else
		error("# infected is 0!")
	end
end

import Base.copy
function copy(stat::Stat)
	return Stat(stat.ns, stat.ne, stat.ni, stat.nr)
end

function add_graph_edges(g::AbstractGraph{T},es::Array{T,2}) where {T<:Unsigned}
	s = size(es)[1]
	@info "Adding graph edges (# edges: ", s, ")"
	oni = Dict{T,T}()
	edges = Array{Tuple{T,T},1}()
	counter = convert(T,1)
	for i in 1:s
		v1 = convert(T,es[i,1])
		v2 = convert(T,es[i,2])
		if !haskey(oni, v1)
			oni[v1] = counter
			counter += convert(T,1)
		end
		if !haskey(oni, v2)
			oni[v2] = counter
			counter += convert(T,1)
		end
		push!(edges, (oni[v1], oni[v2]))
	end

	# add edges
	for edge in edges
		add_edge!(g, edge[1], edge[2])
	end
end

function add_graph_edges(g::AbstractGraph{T},oni::Dict{T,T},es::Array{T,2}) where {T<:Unsigned}
	s = size(es)[1]
	@info "Adding graph edges (# edges: ", s, ")"
	edges = Array{Tuple{T,T},1}()
	for i in 1:s
		v1 = convert(T,es[i,1])
		v2 = convert(T,es[i,2])
		push!(edges, (oni[v1], oni[v2]))
	end
	# add edges
	for edge in edges
		add_edge!(g, edge[1], edge[2])
	end
end

function spread(g::AbstractGraph{T},time::Int,vstats::Dict{T,VState},stat::Stat;e_t::Int=1000,p_i::Float64=0.2,i_t::Int=2000) where {T<:Unsigned}
	@info "Spread infection"
	# E -> I
	for v in vertices(g)
		if vstats[v].status == E::Status && time > vstats[v].time + e_t
			vstats[v].status = I::Status
			vstats[v].time = time
			e2i(stat)
		end
	end
	# I -> R
	for v in vertices(g)
		if vstats[v].status == I::Status && time > vstats[v].time + i_t
			vstats[v].status = R::Status
			vstats[v].time = time
			i2r(stat)
		end
	end
	for edge in edges(g)
		# src infected -> dst exposed
		if vstats[src(edge)].status == I::Status && vstats[dst(edge)].status == S::Status
			if rand() < p_i
				@info("Infecting new node: ", dst(edge))
				# S -> E
				vstats[dst(edge)].status = E::Status
				vstats[dst(edge)].time = time
				s2e(stat)
			end
		end
		# dst infected -> src exposed
		if vstats[dst(edge)].status == I::Status && vstats[src(edge)].status == S::Status
			if rand() < p_i
				@info("Infecting new node: ", src(edge))
				# S -> E
				vstats[src(edge)].status = E::Status
				vstats[src(edge)].time = time
				s2e(stat)
			end
		end
	end
end

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

N_BINS = 4000

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
#vstats[convert(UInt32,ind_pr_max)].status = I::Status
# s2e(stat)

# infect 10 nodes
for j in 1:10
	vstats[convert(UInt32,j)].status = I::Status
	s2e(stat)
	e2i(stat)
end

# set spring layout
l_x, l_y = spring_layout(g_all; C=80)

g_base = SimpleGraph{UInt32}()
add_vertices!(g_base,NV)

# time varies fom 0s to 2418900s to 300s increments (5 min) -> 8063 periods
# 692 vs / 79530 edges
@info("loading bt-detected-internal")
es = Pandas.read_csv("./data/bt-detected-internal.csv")

if SAVE_ANIM
    anim = Animation()
end

## load edge time bins
#ebins = Dict{Int,Array{UInt32,2}}()
#
#@info("loading time bin edges")
#for i in 0:(N_BINS - 1)
#    ess = es[es[:time] < (i+1)*300][es[:time] >= i*300]
#     ea = Array(ess)
#    ebins[i] = ea[1:end .!= 3, 1:end .!= 3]
#end

SA = zeros(Int, N_BINS)
EA = zeros(Int, N_BINS)
IA = zeros(Int, N_BINS)
RA = zeros(Int, N_BINS)

##########
# simulation
##########

for i in 0:(N_BINS - 1)
    # get bin edges
    ess = es[es[:time] < (i+1)*300][es[:time] >= i*300]
    ea = Array(ess)
    eaf =  convert(Array{UInt32,2}, ea[1:end .!= 3, 1:end .!= 3])
    
    # create temporal graph of current interactions
    g = SimpleGraph{UInt32}()
    add_vertices!(g,NV)
    #add_graph_edges(g,oni,ebins[i])
    add_graph_edges(g,oni,eaf)
    #@info "# vs", convert(Int,nv(g))
    #@info "# es", convert(Int,ne(g))
    #@info "----"

    if SAVE_ANIM
    	# enum are 0-based!?
    	# membership = [Int(vstats[v].status)+1 for v in 1:NV]
    	membership = [Int(vstats[v].status) for v in 1:NV]
    	# membership color
    	nodefillc = VCOLOR[membership]

    	p = gplot(g_base, l_x, l_y, nodefillc=nodefillc)
    	output = compose(p)
    
    	j = length(anim.frames) + 1
    	tmpfilename = joinpath(anim.dir, @sprintf("%06d.png", j))
    	Compose.draw(PNG(tmpfilename), output)
    	push!(anim.frames, tmpfilename)
    end
    
    SA[i+1] = stat.ns
    EA[i+1] = stat.ne
    IA[i+1] = stat.ni
    RA[i+1] = stat.nr

    # spread infection
    spread(g,i,vstats,stat,e_t=E_t,p_i=P_i,i_t=I_t)
end

if SAVE_ANIM
	gif(anim, "./plots/animation.gif", fps = 10)
end

Plots.plot(title = "SEIR Simulation", xlabel = "Time", ylabel = "# individuals")
Plots.plot([SA EA IA RA], label = ["S" "E" "I" "R"], linecolor = VCOLOR)

png("./plots/time-series.png")

