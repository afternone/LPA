module LPA

using Graphs
export nsdlpa, nsdlpa1, nmi, voi, modularity, avedegree, similarity, triangle, quadrangle, triquadrsim

include("labelpropagation.jl")
include("modularity.jl")
include("nmi.jl")
include("voi.jl")
include("utils.jl")
include("rankedge.jl")

end # module
