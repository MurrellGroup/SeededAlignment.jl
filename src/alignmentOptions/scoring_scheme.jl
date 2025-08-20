"""
	ScoringScheme(;
		extension_score::Float64=-0.3,
		kmer_length::Int64=18,
		edge_ext_begin::Bool=true,
		edge_ext_end::Bool=true,
		nucleotide_match_score::Float64 = 0.0,
		nucleotide_mismatch_score::Float64 = -0.8,
		nucleotide_score_matrix::Union{Nothing,Matrix{Float64}} = nothing
		codon_match_bonus_score::Float64 = 6.0
)

`ScoringScheme` defines how matches, mismatches, and partially how  gaps are scored during sequence alignment. Beyond
alignment operations governed by `Moveset` it encompasses everything else we can customize about the alignment process. 

This struct is typically passed to functions like `seed_chain_align`, `nw_align` and `msa_codon_align`.

# Fields
-`extension_score::Float64=-0.3`: penalty for affinely extending a gap that is already open. 
-`kmer_length::Int64=18`: length of seeds used in seed_chain_align
-`edge_ext_begin=true`: Whether to subsidize gaps in begining of the alignment
-`edge_ext_end=true`: Whether to subsidize gaps in end of the alignment
-`nucleotide_match_score = 0.0`: score awarded for matching nucleotides (has to be >= 0)
-`nucleotide_mismatch_score = -0.8`: penalty for matching distinct nucleotides (has to be < 0)
-`nucleotide_score_matrix::Union{Nothing, Matrix{Float64}}`: optional custom matrix for scoring nucleotide substitutions.
-`codon_match_bonus_score = 6.0`: score awarded for matching codons

# example
```julia
score_params = ScoringScheme(extension_score = -0.2, mismatch_score = -0.7) # (everything else will be keept at default values)
A = LongDNA{4}("ATGATGATAA")
B = LongDNA{4}("ATGACCCGATAA")
seed_chain_align(A,B scoring=score_params)
```
"""
struct ScoringScheme
	# global extending gap penalty
	extension_score::Float64
	# if seeding - set seed_length
	kmer_length::Int64
	# start and end extension bools
	edge_ext_begin::Bool
	edge_ext_end::Bool
	# substitution scoring:
	nucleotide_match_score::Float64
	nucleotide_mismatch_score::Float64
	# (OPTIONAL) row and column indicies: 1=A, 2=C, 3=G, 4=T
	# Note ScoreScheme is only stack allocated if `nothing` is provided. Providing a matrix can impact performance
	nucleotide_score_matrix::Union{Nothing,Matrix{Float64}}
	# codon_match_bonus_score - only used in ref-query if codon_scoring_on=true in nw_align, seed_chain_align or msa_codon_align
	codon_match_bonus_score::Float64

	# inner constructor 
	function ScoringScheme(extension_score::Float64, kmer_length::Int64, edge_ext_begin::Bool, edge_ext_end::Bool, nucleotide_match_score::Float64, 
		nucleotide_mismatch_score::Float64, nucleotide_score_matrix::Union{Nothing,Matrix{Float64}}, codon_match_bonus_score::Float64)
		# exception handling
		(extension_score < 0) || throw(ArgumentError("extension_score must be negative"))
		(kmer_length > 0) || throw(ArgumentError("kmer_length must be positive integer"))
		# NOTE: Upper and Lower limit might be subject to change
		kmer_length_lower_limit = 6
		kmer_length_upper_limit = 30
		# test kmer_length not too big or too small
		(kmer_length_lower_limit <= kmer_length <= kmer_length_upper_limit) || throw(
			ArgumentError("kmer_length must be within kmer_length limits. Namely we require that:\n
			$(kmer_length_lower_limit) <= kmer_length <= $(kmer_length_upper_limit).")
		)
		# construct ScoringScheme
		new(extension_score, kmer_length, edge_ext_begin, edge_ext_end, nucleotide_match_score, 
			nucleotide_mismatch_score, nucleotide_score_matrix, codon_match_bonus_score)
	end
end

# ScoringScheme construction wrapper
"""
	ScoringScheme(;
		extension_score::Float64=-0.3,
		kmer_length::Int64=18,
		edge_ext_begin::Bool=true,
		edge_ext_end::Bool=true,
		nucleotide_match_score::Float64 = 0.0,
		nucleotide_mismatch_score::Float64 = -0.8,
		nucleotide_score_matrix::Union{Nothing,Matrix{Float64}} = nothing
		codon_match_bonus_score::Float64 = 6.0
)

keyword constructor for `ScoringScheme` with helpful default parameters. It is important to note that
supplying a custom matrix for `nucleotide_score_matrix` might slow performance of some alignment methods.  

# Extended Help

Slow down might occur due to allocating `ScoreScheme` to the heap instead of the stack since
`matrix{Float64}` is dynamically sized. 

# Arguments
-`extension_score::Float64=-0.3`: penalty for affinely extending a gap that is already open. 
-`kmer_length::Int64=18`: length of seeds used in seed_chain_align
-`edge_ext_begin=true`: Whether to subsidize gaps in begining of the alignment
-`edge_ext_end=true`: Whether to subsidize gaps in end of the alignment
-`nucleotide_match_score = 0.0`: score awarded for matching nucleotides (has to be >= 0)
-`nucleotide_mismatch_score = -0.8`: penalty for matching distinct nucleotides (has to be < 0)
-`nucleotide_score_matrix::Union{Nothing, Matrix{Float64}}`: optional custom matrix for scoring nucleotide substitutions.
-`codon_match_bonus_score = 6.0`: score awarded for matching codons

# Returns
-`ScoringScheme`

# example
```julia
score_params = ScoringScheme(extension_score = -0.2, mismatch_score = -0.7) # (everything else will be keept at default values)
A = LongDNA{4}("ATGATGATAA")
B = LongDNA{4}("ATGACCCGATAA")
seed_chain_align(A,B scoring=score_params)

"""
function ScoringScheme(; 
	extension_score::Float64=-0.3, 
	kmer_length::Int64=18, 
	edge_ext_begin=true::Bool, 
	edge_ext_end=true::Bool,
	nucleotide_match_score::Float64 = 0.0,
	nucleotide_mismatch_score::Float64 = -0.8,
	nucleotide_score_matrix=nothing,
	codon_match_bonus_score::Float64 = 6.0,
)
	# call actual constructor
	ScoringScheme(extension_score, kmer_length, edge_ext_begin, edge_ext_end, nucleotide_match_score, 
			nucleotide_mismatch_score, nucleotide_score_matrix, codon_match_bonus_score)
end

import Base: show
function show(io::IO, s::ScoringScheme)
	print(io, "ScoringScheme(",
			  "extension=$(s.extension_score), ",
			  "kmer_length=$(s.kmer_length))",
			  "begin_ext=$(s.edge_ext_begin), ",
			  "end_ext=$(s.edge_ext_end), ",
			  "nucleotide_match_score=$(s.nucleotide_match_score), ",
			  "nucleotide_mismatch_score=$(s.nucleotide_mismatch_score), ",
			  "nucleotide_score_matrix=$(s.nucleotide_score_matrix), ",
			  "codon_match_bonus_score=$(s.codon_match_bonus_score), ")
end

# default ScoringSchemes
"""
	STD_SCORING

constant `ScoreScheme` that stores the default scoring parameters used in alignment methods.

# default parameter values

```julia
const STD_SCORING = ScoringScheme(
	extension_score=-0.3,
	kmer_length=15,
	edge_ext_begin=true,
	edge_ext_end=true,
	nucleotide_mismatch_score = -0.8,
	nucleotide_match_score = 0.0,
	codon_match_bonus_score = 6.0
)

```
"""
const STD_SCORING = ScoringScheme(
	extension_score=-0.3,
	kmer_length=15,
	edge_ext_begin=true,
	edge_ext_end=true,
	nucleotide_mismatch_score = -0.8,
	nucleotide_match_score = 0.0,
	codon_match_bonus_score = 6.0
)