"""
    ScoreScheme

Scoring parameters used for codon-aware sequence alignment.

# Fields
- `match_score::Float64`: Score (typically ≤ 0) awarded for matching codons. Lower is better.
- `mismatch_score::Float64`: Penalty for codon mismatches. Higher values penalize substitutions more strongly.
- `extension_score::Float64`: Cost to extend a gap (indel). Affects how gaps are penalized during alignment.

- `edge_ext_begin::Bool`: If `true`, allows gaps to be extended at the **beginning** of sequences.
- `edge_ext_end::Bool`: If `true`, allows gaps to be extended at the **end** of sequences.

- `kmerlength::Int64`: Length of k-mers used for seeding alignments (if applicable). Ignored if no seeding is used.

# Description
`ScoreScheme` defines how matches, mismatches, and gaps are scored during codon-level sequence alignment. It also controls whether edge gaps can be extended and optionally includes a `kmerlength` for seed-based aligners.

This struct is typically passed to functions like `seed_chain_align` or `msa_codon_align`.
"""
struct ScoreScheme
    match_score::Float64
	mismatch_score::Float64
	extension_score::Float64
	# bools
	edge_ext_begin::Bool
	edge_ext_end::Bool
	# (opt arg) if seeding
	kmerlength::Int64
end


"""
    std_codon_scoring()

Return a standard codon-aware `ScoreScheme` for alignment.

# Returns
- `ScoreScheme`: Default scoring parameters for codon alignment.

# Parameters
- Match score: `0.0`
- Mismatch score: `0.3`
- Gap extension score: `0.1`
- Allow extension at sequence ends: `true` (both ends)
- K-mer length for seeding: `21`

Use this as a default `ScoreScheme` for codon-preserving alignments.
"""
function std_codon_scoring()
	match_score = 0.0
	mismatch_score = 0.3 
	extension_score = 0.1

	edge_ext_begin = true
	edge_ext_end = true

	kmerlength = 21

	return ScoreScheme(match_score,mismatch_score,extension_score,edge_ext_begin,edge_ext_end,kmerlength)
end

function get_all_params(score_params::ScoreScheme)
	return score_params.match_score, score_params.mismatch_score, 
		score_params.extension_score, score_params.edge_ext_begin, score_params.edge_ext_end, score_params.kmerlength
end