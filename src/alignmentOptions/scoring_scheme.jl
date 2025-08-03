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
	# global extending gap penalty
	extension_score::Float64
	# if seeding - set seed_length
	kmer_length::Int64
	# start and end extension bools
	edge_ext_begin::Bool
	edge_ext_end::Bool

	# substitution scoring:
	# row and column indicies: 1=A, 2=C, 3=G, 4=T 
	nucleotide_score_matrix::Matrix{Float64}
	# only used in ref-query if codon_scoring_on=true in nw_align, seed_chain_align or msa_codon_align
	codon_score_matrix::Matrix{Float64}
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
		# TODO add some sort of restriction on scoring matricies. Like matching better than mismatching and not negative
		
		new(extension_score::Float64, kmer_length::Int64, edge_ext_begin::Bool, edge_ext_end::Bool,
			nucleotide_score_matrix::Matrix{Float64}, codon_score_matrix::Matrix{Float64})
	end
end
# TODO definitively need a match/mismatch constructor version because that's easier for people to configure than 20x20 matrix or 4x4

# ScoringScheme construction wrapper # TODO handle stop codon
function ScoringScheme(; extension_score::Float64=-0.4, kmer_length::Int64=18, edge_ext_begin=true::Bool, edge_ext_end=true::Bool, 
		nucleotide_score_matrix=NUC_MATRIX, codon_score_matrix=BLOSUM62)

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
function simple_match_penalty_matrix(match_score, mismatch_score, n=4)
    m = fill(mismatch_score, n, n)
    for i in 1:n
        m[i,i] = match_score
    end
    return m
end

# default movesets and ScoringSchemes

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
# TODO handle stop codon, better score matricies
const STD_SCORING = ScoringScheme(
	extension_score=-0.5,
	kmer_length=15,
	edge_ext_begin=true,
	edge_ext_end=true,
	nucleotide_score_matrix=NUC_MATRIX,
	codon_score_matrix=BLOSUM62
)