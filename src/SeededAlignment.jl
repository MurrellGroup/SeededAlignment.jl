module SeededAlignment

using LinearAlgebra
using BioSequences
using FASTX

# utils
include("utils/fasta_io.jl")
# structs for alignment flexibility
include("alignmentOptions/moveset.jl")
include("alignmentOptions/scoreScheme.jl")
# stats structs for one of the seeding Methods
include("seq_alignment/custom_2D_stats.jl")
# clean_up functionality
include("seq_alignment/clean_frameshifts.jl")
# sequence-sequence alignment
include("seq_alignment/needleman_wunsch.jl")
include("seq_alignment/seed_chain_align.jl")
include("seq_alignment/msa_codon_align.jl")
export 
    # Alignment Methods
    msa_codon_align, seed_chain_align, nw_align,
    # Remove frameshift mutations/noise from a Pairwise Codon Alignment
    clean_frameshifts,
    # Alignment Options types
    Moveset, ScoreScheme, Move
    
end