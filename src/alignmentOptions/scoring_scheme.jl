"""
	ScoringScheme(;
		extension_score::Float64=-0.3,
		kmer_length::Int64=15,
		edge_ext_begin::Bool=true,
		edge_ext_end::Bool=true,
		nucleotide_match_score::Float64 = 0.0,
		nucleotide_mismatch_score::Float64 = -0.8,
		codon_match_bonus_score::Float64 = 6.0
)

....# Description
`ScoringScheme` defines how matches, mismatches, and gaps are scored during nucleotide-level sequence alignment. 
A struct for storing scoring parameters for alignment methods. E.g. nw_align, seed_chain_align and msa_codon_align.

This struct is typically passed to functions like `seed_chain_align`, `nw_align` and `msa_codon_align`.

# Fields
-`extension_score::Float64=-0.3`: 
-`kmer_length::Int64=15`:
-`edge_ext_begin=true`:
-`edge_ext_end=true`:
-`nucleotide_match_score = 0.0`:
-`nucleotide_mismatch_score = -0.8`:
-`codon_match_bonus_score = 6.0`:

# example
```julia
score_params = ScoringScheme(extension_score = -0.2, mismatch_score = -0.7) # (everything else will be keept at default values)
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
	match_nuc_score::Float64
	mismatch_nuc_score::Float64
	# (OPTIONAL) row and column indicies: 1=A, 2=C, 3=G, 4=T
	# Note ScoreScheme is only stack allocated if `nothing` is provided. Providing a matrix can impact performance
	nucleotide_score_matrix::Union{Nothing,Matrix{Float64}} 
	# codon_match_bonus_score - only used in ref-query if codon_scoring_on=true in nw_align, seed_chain_align or msa_codon_align
	codon_match_bonus_score::Float64
	
	codon_score_matrix::Union{Nothing,Matrix{Float64}}
	# codon order: 1=F 2=L 3=S 4=Y 5=C 
	#			   6=W 7=P 8=H 9=Q 10=R 
	#              11=I 12=M 13=T 14=N 15=K
	#              16=V 17=A 18=D 19=E 20=G
	# currently not supported: 21="*", stop codon 22="X" ambigious codon 
	# based on the order of the genetic code

	function ScoringScheme(extension_score::Float64, kmer_length::Int64, edge_ext_begin::Bool, edge_ext_end::Bool,
		nucleotide_score_matrix::Matrix{Float64}, codon_score_matrix::Matrix{Float64})
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
		# TODO force diagional positive on codon_scoring
		# TODO add some sort of restriction on scoring matricies. Like matching better than mismatching and not negative
		
		new(extension_score::Float64, kmer_length::Int64, edge_ext_begin::Bool, edge_ext_end::Bool,
			nucleotide_score_matrix::Matrix{Float64}, codon_score_matrix::Matrix{Float64})
	end
end
# ScoringScheme construction wrapper
function ScoringScheme(; 
	extension_score::Float64=-0.4, 
	kmer_length::Int64=18, 
	edge_ext_begin=true::Bool, 
	edge_ext_end=true::Bool,
	nucleotide_mismatch_score::Float64 = -1.0,
	nucleotide_match_score::Float64 = 0.0,
	codon_mismatch_score::Float64 = 0.0,
	codon_match_score::Float64 = 3.0,
	nucleotide_score_matrix=nothing, 
	codon_score_matrix=nothing
)
	# uses (mis)match_scores if no score matricies are given
	if isnothing(nucleotide_score_matrix)
		nucleotide_score_matrix = simple_match_penalty_matrix(nucleotide_match_score, nucleotide_mismatch_score, 4)
	end
	if isnothing(codon_score_matrix)
		codon_score_matrix = simple_match_penalty_matrix(codon_match_score, codon_match_score, 20)
	end
	# call actual constructor
	ScoringScheme(extension_score::Float64, kmer_length::Int64, edge_ext_begin::Bool, edge_ext_end::Bool,
		nucleotide_score_matrix::Matrix{Float64}, codon_score_matrix::Matrix{Float64})
end

import Base: show
function show(io::IO, s::ScoringScheme)
	print(io, "ScoringScheme(",
			  "extension=$(s.extension_score), ",
			  "begin_ext=$(s.edge_ext_begin), ",
			  "end_ext=$(s.edge_ext_end), ",
			  "kmer_length=$(s.kmer_length))")
end

# matrix constructor helper method
function simple_match_penalty_matrix(match_score::Float64, mismatch_score::Float64, n::Int64)
    m = fill(mismatch_score, n, n)
    for i in 1:n
        m[i,i] = match_score
    end
    return m
end

# default movesets and ScoringSchemes
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
	codon_mismatch_score = 0.0, # TODO this parameter doesn't really work in practice. Works as if always equal 0
	codon_match_score = 6.0
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
	codon_mismatch_score = 0.0, # TODO this parameter doesn't really work in practice. Works as if always equal 0
	codon_match_score = 6.0
)