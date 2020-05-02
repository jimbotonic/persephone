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

using LightGraphs, Pandas
 
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

"""
"""
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

"""
    add_graph_edges(g::AbstractGraph{T},oni::Dict{T,T},es::Array{T,2}) where {T<:Unsigned}

add specified edges to the graph

oni dictionary: old -> new vertex ids

NB: graph g can be directed or undirected
"""
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

"""
    spread(g::AbstractGraph{T},time::Int,vstats::Dict{T,VState},stat::Stat;e_t::Int=1000,p_i::Float64=0.2,i_t::Int=2000) where {T<:Unsigned}
"""
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

"""
    get_bin_edges(::Type{T}, es::DataFrame, time::Int) where {T<:Unsigned}

CSV file headers: s (source), t (target), time

NB: each bin has a 5 min duration (300s)
"""
function get_bin_edges(::Type{T}, es::DataFrame, bin_time::Int) where {T<:Unsigned}
    ess = es[es[:time] < (bin_time+1)*300][es[:time] >= bin_time*300]
    # ess = query(es, :(time >= 0*300 && 1*300 < time))
    #return convert(Array{T,2}, ea[1:end .!= 3, 1:end .!= 3])
    return convert(Array{T,2}, Array(drop(ess, labels=["time"], axis=1)))
end

"""
    get_temporal_graph(::Type{T}, nv::T, oni::Dict{T,T}, edges::Array{T,2}) where {T<:Unsigned}

generate the undirected graph
"""
function get_temporal_graph(::Type{T}, nv::T, oni::Dict{T,T}, edges::Array{T,2}) where {T<:Unsigned}
    g = SimpleGraph{UInt32}()
    add_vertices!(g,nv)
    add_graph_edges(g,oni,edges)
    return g
end

"""
    get_temporal_adj_matrix(nv::Int, es::DataFrame, start::Int, end::Int) 

CSV file headers: source, target, time

NB: each bin has a 5 min duration (300s)
"""
function get_weighted_temporal_adj_matrix(nv::T, es::DataFrame, oni::Dict{T,T}, start_bin::Int, end_bin::Int; is_symmetric=true) where {T<:Unsigned}
    A = zeros(Float64,nv,nv)
    # ess = query(es, :(time >= start_bin*300 && (end_bin+1)*300 < time))
    ess = es[es[:time] < (end_bin+1)*300][es[:time] >= start_bin*300]
    ne = size(ess)[1]
    for i in 1:ne
    	s = convert(T, iloc(ess)[i,1])
    	t = convert(T, iloc(ess)[i,2])
	s = oni[s]
	t = oni[t]
	# increment A entries
	if is_symmetric
        	A[s,t] = A[s,t] + 1.
        	A[t,s] = A[t,s] + 1.
	else
        	A[s,t] = A[s,t] + 1.
	end
    end
    return A
end
