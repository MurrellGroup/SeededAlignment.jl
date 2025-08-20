"""

    seed_chain_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING)

(SeededAlignment wrapper - DE-NOVO)

Computes a heuristically guided global pairwise alignment of two ungapped DNA sequence `A` and `B` based on seeding heuristic. The seeds are then joined together 
by computing an optimal partial alignment between seeds with the Needleman-Wunsch algorithm (nw_align). Optimal in this context meaning with repect to the choosen `Moveset` and `ScoringScheme`.

The advantage of this method is that it is much faster than nw_align and produces similar results for most usecases. 

# Extended Help

# Arguments
- `A::LongDNA{4}`: 1st DNA sequence to be aligned
- `B::LongDNA{4}`: 2nd DNA sequence to be aligned
- `moveset::Moveset=STD_NOISY_MOVESET`: Defines allowable alignment moves (e.g. insertions/deletions and their penalty)
- `scoring::ScoreScheme=STD_SCORING`: Defines alignment scoring together with moveset

# Returns
- `Tuple{LongDNA{4},LongDNA{4}}`: Tuple representation of pairwise alignment of DNA sequences `A` and `B`.

# Example
```julia
# input sequences with no reading frame assumed
A = LongDNA{4}("AATGCTC")
B = LongDNA{4}("ACATGTC")
# produce alignment
alignment = seed_chain_align(A, B)
println(alignment)
#= resulting alignment
alignment = (
	LongDNA{4}("A-ATGCTC"), 
	LongDNA{4}("ACATG-TC")
)
=#
```

"""
function seed_chain_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset = STD_NOISY_MOVESET, scoring::ScoringScheme = STD_SCORING)
    # throw exception if input sequences contains gaps
    !any(x -> x == DNA_Gap, A) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    !any(x -> x == DNA_Gap, B) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    # throw exception if invalid nucleotide letter in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    # check no reference informed moves (Move.ref=true) in moveset
    !contains_ref_move(moveset) || throw(ArgumentError("Moveset contains move(s) that considers reading frame (Move.ref=true)", 
                                                       " when no reference sequence was given!\n","Either set Move.ref=false for all moves in moveset", 
                                                    "if you want to align without a reference sequence. Specify by nw_align(ref=seq1, query=seq2)."))
    # force no clean_up
    do_clean_frameshifts=false
    verbose=false
    # force no codon_scoring_on
    codon_scoring_on=false
    # unpack arguments and call the internal alignment function
    _seed_chain_align(
        A, B, moveset, scoring, codon_scoring_on, do_clean_frameshifts, verbose, false
    )

end

""" 
    seed_chain_align(; 
        	ref::LongDNA{4}, 
        	query::LongDNA{4}, 
        	moveset::Moveset = STD_CODON_MOVESET, 
        	scoring::ScoringScheme = STD_SCORING,
        	codon_scoring_on::Bool = true,
        	do_clean_frameshifts::Bool = false, 
        	verbose::Bool = false)

(SeededAlignment wrapper - CODING given trusted CDS anchor/reference)

Computes a heuristically guided global pairwise alignment of two ungapped CDS (Coding DNA Sequences) `ref` and `query` by using `ref` as an anchor to determine the apprioate reading frame
coordinates for `query` and using a seeding heuristic for speedup. The seeds are then joined together by computing an optimal partial alignment between 
seeds with the Needleman-Wunsch algorithm (nw_align). Optimal in this context means optimal with repect to the choosen `Moveset` and `ScoringScheme`.

The advantage of this method is that it is much faster than nw_align and produces similar results for most usecases. 

# Extended Help

# Arguments
- `ref::LongDNA{4}`: Anchored trusted CDS which decides the reading frame coordinates in the alignment
- `query::LongDNA{4}`: CDS (with possible frameshifts due to e.g. sequencing errors) which is aligned to `ref` and adopts its reading frame coordinates. 
- `moveset::Moveset = STD_CODON_MOVESET`: Defines allowable alignment moves (e.g. insertions/deletions and their penalty)
- `scoring::ScoringScheme = STD_SCORING`: Defines alignment scoring together with moveset
- `codon_scoring_on::Bool = true`: Whether to apply additional scoring on codon-level 
- `do_clean_frameshifts::Bool = false`: Whether to clean the alignment output of gaps which cause frameshifts - this produces a protein alignment on a nucleotide level. 
- `verbose::Bool = false`: Whether to verbosely display what edits were made during the cleaning of frameshifts. 

# Returns
- `Tuple{LongDNA{4},LongDNA{4}}`: Tuple representation of pairwise alignment of DNA sequences `ref` and `query`. 
Note that this represents a protein alignment on a nucleotide level if `(do_clean_frameshifts == true)`.

# Example
```julia
anchor_CDS =    LongDNA{4}("ATGCCAGTA")
# untrusted_CDS may contain some frameshift errors due to e.g. sequencing or annotation errors.
untrusted_CDS = LongDNA{4}("ATGTA") 
# frameshift errors are removed from the cleaned alignment
cleaned_CDS_alignment = seed_chain_align(ref=anchor_CDS, query=untrusted_CDS, clean_frameshifts=true)
println(cleaned_CDS_alignment)
#= resulting alignment:

cleaned_CDS_alignment = (
	LongDNA{4}("ATGCCAGTA"), 
	LongDNA{4}("ATG---NTA")
)
Here 'N' denotes ambigious nucleotide.
=#
```

"""
function seed_chain_align(; 
    ref::LongDNA{4},
    query::LongDNA{4},
    moveset::Moveset = STD_CODON_MOVESET,
    scoring::ScoringScheme = STD_SCORING,
    codon_scoring_on::Bool = true,
    do_clean_frameshifts::Bool = false,
    verbose::Bool = false)
    # throw exception if input sequences contains gaps
    !any(x -> x == DNA_Gap, ref) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    !any(x -> x == DNA_Gap, query) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    # throw exception if invalid nucleotide letter in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), ref) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), query) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    # check that moveset takes reference reading frame into account
    contains_ref_move(moveset) || throw(ArgumentError("Invalid Moveset for reference to query alignment!\n
                                        At least one Move in Moveset must consider reference reading (Move.ref=true)
                                         - in other words codon insertions or deletions must be allowed."))
    # unpack arguments and call the internal alignment function
    _seed_chain_align(
        ref, query, moveset, scoring, codon_scoring_on, do_clean_frameshifts, verbose, true
    )
end

@inbounds function _seed_chain_align(A::LongDNA{4}, B::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING, 
    codon_scoring_on=true::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool, A_is_ref::Bool=false) #; seed_debug_mode=true)

    # unpack parameters
    k = scoring.kmer_length
    vgap_moves = moveset.vert_moves
    hgap_moves = moveset.hor_moves
    nucleotide_score_matrix = scoring.nucleotide_score_matrix
    match_score = scoring.nucleotide_match_score
    mismatch_score = scoring.nucleotide_mismatch_score
    codon_match_bonus_score = scoring.codon_match_bonus_score
    extension_score = scoring.extension_score
    edge_ext_begin = scoring.edge_ext_begin
    edge_ext_end = scoring.edge_ext_end
    # seeding heuristic
    kmerMatches = find_kmer_matches(A, B, k, A_is_ref=A_is_ref)
    # two possible seeding heuristics
    kmerPath = select_kmer_path(kmerMatches, length(A), length(B), nucleotide_score_matrix, match_score, mismatch_score, 
        vgap_moves, hgap_moves, extension_score, k)
    #kmerPath = select_max_correlation_kmer_path(kmerMatches, k)
    # parameters for joining kmers/seeds
    extra_kmer_margin = 6
    # Join kmers using needleman-wunsch
    k = k-extra_kmer_margin
    prevA = -k+1
    prevB = -k+1
    result = [LongDNA{4}(""), LongDNA{4}("")]
    @views for kmer in kmerPath
        # TODO we want to do "result .*=" all in one step preferably
        if !(kmer.posA == prevA + k && kmer.posB == prevB + k)
            if prevA == -k+1 && prevB == -k+1
                # align without cleaning frameshifts
                alignment = _nw_align(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], vgap_moves, hgap_moves, extension_score, 
                    edge_ext_begin, false, match_score, mismatch_score, nucleotide_score_matrix, codon_match_bonus_score, codon_scoring_on)
                # add alignment unless seed_debug_mode
                #= (only for internal debug)
                if seed_debug_mode
                    obfuscate_nucleotides!.(alignment)
                end
                =#
                result .*= alignment
            else
                # align without cleaning frameshifts
                alignment = _nw_align(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], vgap_moves, hgap_moves, extension_score, 
                    false, false, match_score, mismatch_score, nucleotide_score_matrix, codon_match_bonus_score, codon_scoring_on)
                # add alignment unless seed_debug_mode
                #= (only for internal debug) 
                if seed_debug_mode
                    obfuscate_nucleotides!.(alignment)
                end
                =#
                result .*= alignment 
            end
        end
        result .*= [A[kmer.posA : kmer.posA + k - 1], B[kmer.posB : kmer.posB + k - 1]]
        prevA = kmer.posA
        prevB = kmer.posB
    end
    # align the remaining parts of the sequences - without cleaning frameshifts
    alignment = _nw_align(A[prevA + k : end], B[prevB + k : end], vgap_moves, hgap_moves, extension_score, 
        false, edge_ext_end, match_score, mismatch_score, nucleotide_score_matrix, codon_match_bonus_score, codon_scoring_on)
    # add alignment unless seed_debug_mode
    #= (only for internal debug) 
    if seed_debug_mode
        obfuscate_nucleotides!.(alignment)
    end
    =#
    result .*= alignment 
    # clean_frameshifts 
    if do_clean_frameshifts # && !seed_debug_mode (only for internal debug)
        result[1], result[2] = clean_frameshifts(result[1], result[2], verbose = verbose)
    end
    # return result as Tuple
    return result[1], result[2]
end

#= (only for internal debug)
function obfuscate_nucleotides!(seq::LongDNA{4})
    for i in eachindex(seq)
        if seq[i] != DNA_Gap && seq[i] != LongDNA{4}("N")
            seq[i] = DNA_R        
        end
    end
end
=#