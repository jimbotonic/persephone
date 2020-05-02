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

using LightGraphs, Clustering

using Adjacently.io
using Adjacently.graph

include("./utils.jl")

# Documentation
# https://juliastats.org/Clustering.jl/stable/hclust.html

##########
# simulation variables
##########

# whole period Copenhagen graph
NV = convert(UInt32, 692)
NE = 79530


@info("loading bt-detected-internal-notime")
g_all = SimpleGraph{UInt32}()
oni = load_adjacency_list_from_csv(UInt32, g_all, "./data/bt-detected-internal-notime.csv")
# write_mgs3_graph(g_all, "bt-detected-internal-notime")

A = get_sparse_adj_matrix(g_all) 
mres = mcl(A)

# 3 clusters
@info("MCL A: # clusters ", length(mres.counts))
@info("Counts: ", mres.counts)
@info("Assignements: ", mres.assignments)

# time varies fom 0s to 2418900s to 300s increments (5 min) -> 8063 periods
# 692 vs / 79530 edges
# CSV file headers: source, target, time
@info("loading bt-detected-internal")
es = Pandas.read_csv("./data/bt-detected-internal.csv")

A2 = get_weighted_temporal_adj_matrix(NV, es, oni, 0, 100) 
mres2 = mcl(A2)

@info("MCL A2: # clusters ", length(mres2.counts))
@info("Counts: ", mres2.counts)
@info("Assignements: ", mres2.assignments)

hres2 = hclust(A2)
ha = cutree(hres2, k=2)

@info("H assignments: ", ha)


