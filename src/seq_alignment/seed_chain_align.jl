"""
    seed_chain_align(A::LongDNA{4},B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING
    
seed_chain_align wrapper - no reference, i.e. makes no assumptions about the two sequences. 

Computes a heuristically guided global pairwise alignment of two ungapped `LongDNA{4}` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to the `Moveset` and `ScoringScheme`.
"""
function seed_chain_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING)
    # check no reference informed moves (Moev.ref=true) in moveset
    !contains_ref_move(moveset) || throw(ArgumentError("Moveset contains move(s) that considers reading frame (Move.ref=true)", 
                                                       " when no reference sequence was given!\n","Either set Move.ref=false for all moves in moveset", 
                                                    "if you want to align without a reference sequence. Specify by nw_align(ref=seq1, query=seq2)."))
    # force no clean_up
    do_clean_frameshifts=false
    verbose=false
    # force no match_codons
    match_codons=false
    # unpack arguments and call the internal alignment function
    seed_chain_align(
        A, B, moveset.vert_moves, moveset.hor_moves,
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_score_matrix, 
        scoring.kmer_length, match_codons, do_clean_frameshifts, verbose
    )

end

""" 
    
    seed_chain_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING,
        do_clean_frameshifts=false::Bool, verbose=false::Bool, match_codons=true::Bool)

seed_chain_align wrapper - reference informed, i.e. assumes one of the sequence has intact reading frame. 

Heuristically aligns a `query` sequence to a `ref` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to a codon-aware `Moveset` and `ScoringScheme`.

**NOTE** We always assume the readingFrame is 1
"""
function seed_chain_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING, 
    match_codons=true::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool)

    # check that moveset takes reference reading frame into account
    contains_ref_move(moveset) || throw(ArgumentError("Invalid Moveset for reference to query alignment!\n", 
                                        "At least one Move in Moveset must consider reference reading (Move.ref=true)",
                                        " - in other words codon insertions or deletions must be allowed."))
    # unpack arguments and call the internal alignment function
    seed_chain_align(
        ref, query, moveset.vert_moves, moveset.hor_moves,
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_score_matrix,
        scoring.kmer_length, match_codons, do_clean_frameshifts, verbose
    )
end

function seed_chain_align(A::LongDNA{4}, B::LongDNA{4}, vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, 
    match_score_matrix::Matrix{Float64}, extension_score::Float64, codon_score_matrix::Matrix{Float64}=BLOSUM62, kmer_length::Int64 = 18,
    match_codons=false::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool) where {X, Y}

    # throw exception if invalid alphabet in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    # Abbreviations
    k = kmer_length
    # seeding heuristic
    kmerMatches = find_kmer_matches(A, B, k)
    kmerPath = select_kmer_path(kmerMatches, length(A), length(B), match_score_matrix, vgap_moves, hgap_moves, extension_score, k)    
    # parameters for joining kmers/seeds
    extra_kmer_margin = 6
    # Join kmers using needleman-wunsch
    k = k-extra_kmer_margin
    prevA = -k+1
    prevB = -k+1
    result = [LongDNA{4}(""), LongDNA{4}("")]
    for kmer in kmerPath
        # shift start of kmer to start of next codon
        offset_from_codon_boundary = ((kmer.posA) % 3)
        kmer.posA += (3-offset_from_codon_boundary)+4
        kmer.posB += (3-offset_from_codon_boundary)+4
        if !(kmer.posA == prevA + k && kmer.posB == prevB + k)
            if prevA == -k+1 && prevB == -k+1
                # align without cleaning frameshifts
                alignment = nw_align(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], vgap_moves, hgap_moves, 
                    match_score_matrix, extension_score, codon_score_matrix, true, false, match_codons
                )
                result .*= alignment
            else
                # align without cleaning frameshifts
                alignment = nw_align(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], vgap_moves, hgap_moves, 
                    match_score_matrix, extension_score, codon_score_matrix, false, false, match_codons
                )
                result .*= alignment
            end
        end
        result .*= [A[kmer.posA : kmer.posA + k - 1], B[kmer.posB : kmer.posB + k - 1]]
        prevA = kmer.posA
        prevB = kmer.posB
    end
    # align the remaining parts of the sequences - without cleaning frameshifts
    result .*= nw_align(A[prevA + k : end], B[prevB + k : end], vgap_moves, hgap_moves, 
        match_score_matrix, extension_score, codon_score_matrix, false, true, match_codons
    )
    if do_clean_frameshifts
        result[1], result[2] = clean_frameshifts(result[1], result[2], verbose = verbose)
    end
    # return result as Tuple
    return result[1], result[2]
end