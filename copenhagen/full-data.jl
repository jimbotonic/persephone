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

#using Pkg
#Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs, GraphPlot, Compose, Cairo, Fontconfig

using Adjacently.util
using Adjacently.io
using Adjacently.graph

function load_dataset(input_path::AbstractString,output_filename::AbstractString)
	g = SimpleGraph{UInt32}()
	load_adjacency_list_from_csv(UInt32, g, input_path)
	
	@info("# vertices:", convert(Int,nv(g)))
	@info("# edges:", ne(g))
	
	write_mgs3_graph(g, output_filename)
	# write_mgs4_graph(g, output_filename)

	layout=(args...)->spring_layout(args...; C=40)
	draw(PNG("g.png", 16cm, 16cm), gplot(g, layout=layout))
end

@info("loading bt-detected-internal-notime")
load_dataset("./data/bt-detected-internal-notime.csv", "bt-detected-internal-notime")
