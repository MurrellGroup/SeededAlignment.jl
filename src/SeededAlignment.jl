module SeededAlignment

using LinearAlgebra
using BioSequences

export Move, nw_align, seed_chain_align

include("needleman_wunsch.jl")
include("custom_2D_stats.jl")
include("seed_chain_align.jl")

end
