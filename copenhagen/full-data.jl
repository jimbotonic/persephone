using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

include("../src/util.jl")
include("../src/io.jl")
include("../src/graph.jl")

using GraphPlot, Compose, Cairo, Fontconfig

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
