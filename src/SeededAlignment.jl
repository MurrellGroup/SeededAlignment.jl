module SeededAlignment


using FASTX
# export dna sequence datatype
using BioSequences
using BioSequences: LongDNA
export LongDNA
# exports utilities
include("utils/fasta_io.jl")
# structs for alignment flexibility
include("alignmentOptions/constants_and_helpers.jl")
include("alignmentOptions/moveset.jl")
include("alignmentOptions/scoring_scheme.jl")
# clean_up functionality for pairwise and multiple sequence alignment
include("seq_alignment/clean_frameshifts.jl")
# sequence-sequence alignment
include("seq_alignment/needleman_wunsch.jl")
include("seq_alignment/seeding.jl")
include("seq_alignment/seed_chain_align.jl")
include("seq_alignment/msa_codon_align.jl")
export 
    # Alignment Methods
    msa_codon_align, seed_chain_align, nw_align,
    # Remove frameshift mutations/noise from a Pairwise Codon Alignment
    clean_frameshifts,
    # Alignment Options types
    Moveset, ScoringScheme, Move,
    # utilities
    read_fasta, write_fasta
end