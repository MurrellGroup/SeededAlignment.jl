"""
```julia
ScoringScheme(; match_score=0.0, mismatch_score=0.5,extension_score=0.1,edge_ext_begin=true,edge_ext_end=true,kmerlength=21)
```

A struct for storing scoring parameters used for sequence alignment. Uses a keyword constructor. 

# Fields
- `match_score::Float64`: Score (typically ≤ 0) awarded for matching nucleotide. Lower is better.
- `mismatch_score::Float64`: Penalty for nucleotide mismatches. Higher values penalize substitutions more strongly.
- `extension_score::Float64`: Cost to extend a gap (indel). Affects how gaps are penalized during alignment.

- `edge_ext_begin::Bool`: If `true`, allows gaps to be extended at the **beginning** of sequences.
- `edge_ext_end::Bool`: If `true`, allows gaps to be extended at the **end** of sequences.

- `kmerlength::Int64`: Length of k-mers used for seeding alignments (if applicable). Ignored if no seeding is used.

# Description
`ScoringScheme` defines how matches, mismatches, and gaps are scored during nucleotide-level sequence alignment. 

This struct is typically passed to functions like `seed_chain_align` or `msa_codon_align`.

# example
```julia
score_params = ScoringScheme(extension_score = 0.3, mismatch_score = 0.7) # (everything else will be keept at default values)
```
"""
struct ScoringScheme
	# 1=A, 2=C, 3=G, 4=T for both row and column indicies
	nucleotide_score_matrix::Matrix{Float64}
	# extending gap penalty
	extension_score::Float64
	# start and end extension bools 
	edge_ext_begin::Bool
	edge_ext_end::Bool
	# (opt arg) if match_codons in aligner function
	# TODO replace with codon scoring matrix
	codon_match_bonus::Float64
	# (opt arg) if seeding
	kmerlength::Int64
	
end

# ScoringScheme construction wrapper
function ScoringScheme(; nucleotide_score_matrix::Union{Matrix{Float64},Nothing}=nothing, match_score::Float64 = 0.0, mismatch_score::Float64=0.7,
		extension_score::Float64=0.4, edge_ext_begin=true::Bool, edge_ext_end=true::Bool, 
		codon_match_bonus::Float64=-2.0, kmerlength::Int64=18)

	# give a default simple match/mismatch matrix if no score_matrix provided
	if nucleotide_score_matrix == nothing
		nucleotide_score_matrix = simple_match_penalty_matrix(match_score, mismatch_score)
	end
	# call actual constructor
	return ScoringScheme(nucleotide_score_matrix,extension_score,
		edge_ext_begin,edge_ext_end,
		codon_match_bonus,kmerlength
	)	
end

import Base: show
function show(io::IO, s::ScoringScheme)
	print(io, "ScoringScheme(",
			  "match=$(s.match_score), ",
			  "mismatch=$(s.mismatch_score), ",
			  "extension=$(s.extension_score), ",
			  "begin_ext=$(s.edge_ext_begin), ",
			  "end_ext=$(s.edge_ext_end), ",
			  "codon_match_bonus=$(s.codon_match_bonus), ",
			  "kmer=$(s.kmerlength))")
end

"""
    std_codon_scoring()

Return a standard codon-aware `ScoringScheme` for alignment.

# Returns
- `ScoringScheme`: Default scoring parameters for codon alignment.

# Parameters
- Match score: `0.0`
- Mismatch score: `0.3`
- Gap extension score: `0.1`
- Allow extension at sequence ends: `true` (both ends)
- K-mer length for seeding: `18`

Use this as a default `ScoringScheme` for codon-preserving alignments.
"""
function std_scoring()
	# params for nucleotide_score_matrix
	match_score = 0.0
	mismatch_score = 0.7
	# scoring stuff
	extension_score = 0.4
	# start and end extension bools
	edge_ext_begin = true
	edge_ext_end = true
	# (opt arg) if match_codons in aligner function
	codon_match_bonus = -2.0
	# (opt arg) if seeding
	kmerlength = 18

	return ScoringScheme(match_score=match_score,
					   mismatch_score=mismatch_score,
					   extension_score=extension_score,
					   edge_ext_begin=edge_ext_begin,
					   edge_ext_end=edge_ext_end,
					   codon_match_bonus=codon_match_bonus,
					   kmerlength=kmerlength
	)
end