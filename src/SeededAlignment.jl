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
include("seq_alignment/clean_alignment_readingframe.jl")
# sequence-sequence alignment
include("seq_alignment/needleman_wunsch.jl")
include("seq_alignment/seed_chain_align.jl")
include("seq_alignment/msa_codon_align.jl")
export 
    # Alignment Methods
    msa_codon_align, seed_chain_align, nw_align,
    # Alignment Options constructors
    MoveSet, ScoreScheme, Move, 
    # Remove single indel noise from a Pairwise Codon Alignment
    clean_alignment_readingframe,
    # Alignment Options Helper functions
    std_codon_moveset, std_codon_scoring, pairwise_noisy_moveset
end