"""
    seed_chain_align(A::LongDNA{4},B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING
    
seed_chain_align wrapper - no reference, i.e. makes no assumptions about the two sequences. 

Computes a heuristically guided global pairwise alignment of two ungapped `LongDNA{4}` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to the `Moveset` and `ScoringScheme`.
"""
function seed_chain_align(A::LongDNA{4}, B::LongDNA{4};
    moveset::Moveset = STD_NOISY_MOVESET,
    scoring::ScoringScheme = STD_SCORING
)
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
    seed_debug_mode=false
    _seed_chain_align(
        A, B, moveset, scoring, codon_scoring_on, do_clean_frameshifts, verbose, seed_debug_mode=seed_debug_mode
    )

end

""" 
    
    seed_chain_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING,
        do_clean_frameshifts=false::Bool, verbose=false::Bool, codon_scoring_on=true::Bool)

seed_chain_align wrapper - reference informed, i.e. assumes one of the sequence has intact reading frame. 

Heuristically aligns a `query` sequence to a `ref` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to a codon-aware `Moveset` and `ScoringScheme`.
"""
function seed_chain_align(; 
    ref::LongDNA{4},
    query::LongDNA{4},
    moveset::Moveset = STD_CODON_MOVESET,
    scoring::ScoringScheme = STD_SCORING,
    codon_scoring_on::Bool = true,
    do_clean_frameshifts::Bool = false,
    verbose::Bool = false
)
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
    seed_debug_mode=false
    _seed_chain_align(
        ref, query, moveset, scoring, codon_scoring_on, do_clean_frameshifts, verbose, seed_debug_mode=seed_debug_mode
    )
end

@inbounds function _seed_chain_align(A::LongDNA{4}, B::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING, 
    codon_scoring_on=true::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool; seed_debug_mode=false::Bool)

    # unpack parameters
    k = scoring.kmer_length
    vgap_moves = moveset.vert_moves
    hgap_moves = moveset.hor_moves
    match_score_matrix = scoring.nucleotide_score_matrix
    extension_score = scoring.extension_score
    codon_score_matrix = scoring.codon_score_matrix
    edge_ext_begin = scoring.edge_ext_begin
    edge_ext_end = scoring.edge_ext_end
    # seeding heuristic
    kmerMatches = find_kmer_matches(A, B, k)
    # two possible seeding heuristics # TODO compare performance of these
    kmerPath = select_kmer_path(kmerMatches, length(A), length(B), match_score_matrix, vgap_moves, hgap_moves, extension_score, k)
    #kmerPath = select_max_correlation_kmer_path(kmerMatches, k)
    # parameters for joining kmers/seeds
    extra_kmer_margin = 6
    # Join kmers using needleman-wunsch
    k = k-extra_kmer_margin
    prevA = -k+1
    prevB = -k+1
    result = [LongDNA{4}(""), LongDNA{4}("")]
    @views for kmer in kmerPath
        # TODO make so seeding considers codons for better performance, especially then we can combine kmers better
        # TODO we want to do "result .*=" all in one step preferably
        # shift start of kmer to start of next codon
        offset_from_codon_boundary = ((kmer.posA) % 3)
        new_posA = (3-offset_from_codon_boundary)+4 + kmer.posA
        new_posB = (3-offset_from_codon_boundary)+4 + kmer.posB
        if !(new_posA == prevA + k && new_posB == prevB + k)
            if prevA == -k+1 && prevB == -k+1
                # align without cleaning frameshifts
                alignment = _nw_align(A[prevA + k : new_posA - 1], B[prevB + k : new_posB - 1], vgap_moves, hgap_moves, 
                    match_score_matrix, extension_score, codon_score_matrix, edge_ext_begin, false, codon_scoring_on
                )
                # add alignment unless seed_debug_mode
                if seed_debug_mode
                    obfuscate_nucleotides!.(alignment)
                end
                result .*= alignment
            else
                # align without cleaning frameshifts
                alignment = _nw_align(A[prevA + k : new_posA - 1], B[prevB + k : new_posB - 1], vgap_moves, hgap_moves, 
                    match_score_matrix, extension_score, codon_score_matrix, false, false, codon_scoring_on
                )
                # add alignment unless seed_debug_mode
                if seed_debug_mode
                    obfuscate_nucleotides!.(alignment)
                end
                result .*= alignment 
            end
        end
        result .*= [A[new_posA : new_posA + k - 1], B[new_posB : new_posB + k - 1]]
        prevA = new_posA
        prevB = new_posB
    end
    # align the remaining parts of the sequences - without cleaning frameshifts
    alignment = _nw_align(A[prevA + k : end], B[prevB + k : end], vgap_moves, hgap_moves, 
        match_score_matrix, extension_score, codon_score_matrix, false, edge_ext_end, codon_scoring_on
    )
    # add alignment unless seed_debug_mode
    if seed_debug_mode
        obfuscate_nucleotides!.(alignment)
    end
    result .*= alignment 
    # clean_frameshifts 
    if do_clean_frameshifts && !seed_debug_mode
        result[1], result[2] = clean_frameshifts(result[1], result[2], verbose = verbose)
    end
    # return result as Tuple
    return result[1], result[2]
end
# only for internal debug
function obfuscate_nucleotides!(seq::LongDNA{4})
    for i in eachindex(seq)
        if seq[i] != DNA_Gap && seq[i] != LongDNA{4}("N")
            seq[i] = DNA_R        
        end
    end
end