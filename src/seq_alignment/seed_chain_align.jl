"""
    seed_chain_align(A::LongDNA{4},B::LongDNA{4}; moveset::Moveset=std_noisy_moveset(), scoring::ScoringScheme=std_scoring()
    
seed_chain_align wrapper - no reference, i.e. makes no assumptions about the two sequences. 

Computes a heuristically guided global pairwise alignment of two ungapped `LongDNA{4}` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to the `Moveset` and `ScoringScheme`.
"""
function seed_chain_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=std_noisy_moveset(), scoring::ScoringScheme=std_scoring())
    # force no clean_up
    do_clean_frameshifts=false
    verbose=false
    # force no match_codons
    match_codons=false
    # unpack arguments and call the internal alignment function
    seed_chain_align(
        A, B, moveset.vert_moves, moveset.hor_moves,
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_match_bonus, 
        scoring.kmerlength, match_codons, do_clean_frameshifts, verbose
    )

end

""" 
    
    seed_chain_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=std_codon_moveset(), scoring::ScoringScheme=std_scoring(),
        do_clean_frameshifts=false::Bool, verbose=false::Bool, match_codons=true::Bool)

seed_chain_align wrapper - reference informed, i.e. assumes one of the sequence has intact reading frame. 

Heuristically aligns a `query` sequence to a `ref` sequence based on seeding. The seeds are then joined together 
by doing a partial alignment with nw_align between seeds optimally with repect to a codon-aware `Moveset` and `ScoringScheme`.

**NOTE** We always assume the readingFrame is 1
"""
function seed_chain_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=std_codon_moveset(), scoring::ScoringScheme=std_scoring(), 
    match_codons=true::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool)
    # unpack arguments and call the internal alignment function
    seed_chain_align(
        ref, query, moveset.vert_moves, moveset.hor_moves,
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_match_bonus,
        scoring.kmerlength, match_codons, do_clean_frameshifts, verbose
    )
end

function find_kmer_matches(A::LongDNA{4}, B::LongDNA{4}, kmerLength)

    # Abbreviations
    k = kmerLength
    m = length(A)
    n = length(B)

    # List of all kmer matches between A and B
    kmerMatches = KmerMatch[]
    # Produce dictionary of all kmers in A
    kmerDict = Dict{LongDNA{4}, Array{Int64}}()
    for i in 1:m-k+1
        kmer = A[i:i+k-1]
        if kmer in keys(kmerDict)
            push!(kmerDict[kmer], i)
        else
            kmerDict[kmer] = [i]
        end
    end

    # Search B for any matching kmers
    diagonals = fill(typemin(Int64), m+n) # diagonals[i] = the rightmost kmer in diagonal i. (used to find overlaps in matches)
    repetitionTolerance = 5 # Threshold for ignoring repeted kmers
    for iB in 1 : n-k+1
        kmer = B[iB : iB+k-1]
        if (haskey(kmerDict, kmer) && length(kmerDict[kmer]) <=  repetitionTolerance)
            
            #Add found match(es) to list
            for iA in kmerDict[kmer]
                if diagonals[iA - iB + n + 1] + k <= iA
                    push!(kmerMatches, KmerMatch(iA, iB))
                    diagonals[iA - iB + n + 1] = iA
                end
            end
        end
    end

    return kmerMatches
end

function select_kmer_path(kmerMatches, m::Int64, n::Int64, match_score_matrix::Matrix{Float64}, 
    vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, extension_score::Float64, k::Int64) where {X, Y}
    
    # Produce constants used for estimating scores without A and B
    min_match_score = minimum(t -> match_score_matrix[t, t], 1 : 4)
    gap_score_estimate = minimum(move -> move.score / move.step, (vgap_moves..., hgap_moves...))
    if gap_score_estimate > extension_score >= 0
        gap_score_estimate = (gap_score_estimate + 2 * extension_score) / 3 # approx_gap_score undefined, replaced with gap_score_estimate
    end
    gap_score_estimate -= min_match_score/2
    match_score_estimate = sum(match_score_matrix) / 16 - min_match_score

    matchCount = length(kmerMatches)
    
    #Produce list of endpoints - two for each kmer
    endpoints = Array{Endpoint}(undef,matchCount * 2)
    for i in 1 : matchCount
        match = kmerMatches[i] 
        endpoints[2 * i - 1] = Endpoint(match.posA, match.posB, i, true)
        endpoints[2 * i] = Endpoint(match.posA + k, match.posB + k, i, false)
    end
    sort!(endpoints, by = e -> (e.x, e.y, e.isBeginning))
    
    # Transitions
    bestConnections = Array{Int64}(undef, matchCount)
    bestScores = Array{Float64}(undef, matchCount)

    #Chains speed up the DP
    chainNum = 0
    chains = Vector{Vector{Int64}}() # Disjoint sequences of compatible kmers
    chainIds = Array{Int64}(undef, matchCount) # The chain of each kmer
    occupiedChains = Vector{Bool}() # chains with recently added kmers
    chainLinks = Vector{Vector{Int64}}() # chainLinks[x,y] = the last kmer in chain y which some kmer in chain x can connect with.
    
    #Find best path through kmers
    for e in endpoints
        if e.isBeginning #Begin kmer
            
            #Find the closest chain it can connect with
            bestChain = -1
            bestChainY = -1
            for j in 1 : chainNum
                if !occupiedChains[j]
                    chainY = kmerMatches[last(chains[j])].posB + k
                    if chainY < e.y && chainY > bestChainY
                        bestChain = j
                        bestChainY = chainY
                    end
                end
            end

            #Add it to the chain or create a new one
            if bestChain == -1
                chainNum += 1
                bestChain = chainNum
                push!(chainLinks, fill(0, chainNum)) 
                for j in 1 : chainNum-1
                    push!(chainLinks[j], 0)
                end
                push!(chains, [])
                push!(occupiedChains, true)
            else
                occupiedChains[bestChain] = true
            end
            
            #Find best connection
            chainIds[e.id] = bestChain
            bestScores[e.id] = getConnectionScore(KmerMatch(-k,-k), KmerMatch(e.x, e.y), k, gap_score_estimate, match_score_estimate)
            bestConnections[e.id] = e.id
            for j in 1 : chainNum

                # Update chain links
                t = chainLinks[bestChain][j]
                while t + 1 <= length(chains[j]) && kmerMatches[chains[j][t + 1]].posB + k <= e.y
                    t += 1
                end
                chainLinks[bestChain][j] = t

                # Pick best chain link
                if t != 0
                    candidate = chains[j][t]
                    score = getConnectionScore(kmerMatches[candidate], KmerMatch(e.x, e.y), k, gap_score_estimate, match_score_estimate) + bestScores[candidate]
                    if score < bestScores[e.id]
                        bestScores[e.id] = score
                        bestConnections[e.id] = candidate
                    end
                end
            end
        else
            #End kmer
            occupiedChains[chainIds[e.id]] = false
            push!(chains[chainIds[e.id]], e.id)
        end
    end
    
    #Back propagation
    kmerPath = Vector{KmerMatch}() #Final kmer choice
    if matchCount > 0
        kmer = argmin([bestScores[i] + getConnectionScore(kmerMatches[i], KmerMatch(m, n), k, gap_score_estimate, match_score_estimate) for i in 1 : matchCount])
        while true
            push!(kmerPath, kmerMatches[kmer])
            if bestConnections[kmer] == kmer
                break
            else
                kmer = bestConnections[kmer]
            end
        end
        reverse!(kmerPath)
    end

    return kmerPath
end

function seed_chain_align(A::LongDNA{4}, B::LongDNA{4}, vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, 
    match_score_matrix::Matrix{Float64}, extension_score::Float64 = -1.0, codon_match_bonus::Float64 = -2.0, kmerLength::Int64 = 12,
    match_codons=false::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool) where {X, Y}

    # throw exception if invalid alphabet in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    # Abbreviations
    k = kmerLength
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
                    match_score_matrix, extension_score, codon_match_bonus, true, false, match_codons
                )
                result .*= alignment
            else
                # align without cleaning frameshifts
                alignment = nw_align(A[prevA + k : kmer.posA - 1], B[prevB + k : kmer.posB - 1], vgap_moves, hgap_moves, 
                    match_score_matrix, extension_score, codon_match_bonus, false, false, match_codons
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
        match_score_matrix, extension_score, codon_match_bonus, false, true, match_codons
    )
    if do_clean_frameshifts
        result[1], result[2] = clean_frameshifts(result[1], result[2], verbose = verbose)
    end
    # return result as Tuple
    return result[1], result[2]
end

# helper types and functions
function CorrelationScore(s::SampleMetrics)
    return sqrt(s.n) * s.correlation
end

mutable struct KmerMatch
    posA::Int64
    posB::Int64
end

# Estimate pealty between kmer matches x and y
function getConnectionScore(x::KmerMatch, y::KmerMatch, k::Int64, gap_score_estimate::Float64, match_score_estimate::Float64)
    m = y.posA - x.posA - k
    n = y.posB - x.posB - k
    return abs(m - n) * gap_score_estimate + min(m, n) * match_score_estimate
end

struct Endpoint
    x::Int64
    y::Int64
    id::Int64
    isBeginning::Bool
end

#Kmer selection variant algorithm (unused currently) will be compared in benchmark
function select_max_correlation_kmer_path(kmerMatches)

    # Initialize
    matchCount = length(kmerMatches)
    remaining_matches = matchCount
    deletionFlags = BitArray([0 for i in 1:matchCount])
    kmerMetrics = SampleMetrics(map(x -> x.posA, kmerMatches), map(x -> x.posB, kmerMatches))
    corScore = CorrelationScore(kmerMetrics)
    deletion_candidates = Array{Tuple{Float64, Int64}}(undef, matchCount)
    pruning_steps = 10
    max_pruning_fraction = 0.15
    
    # Prune kmers
    for i in 1 : pruning_steps
        if kmerMetrics.n <= 2
            break
        end
        for j in 1 : matchCount
            if !deletionFlags[j]
                kmer = kmerMatches[j]
                deletion_candidates[j] = (CorrelationScore(remove_point(kmerMetrics, kmer.posA, kmer.posB)), j)
            end
        end
        max_pruning_num = floor(Int64,remaining_matches*max_pruning_fraction + 1.0)
        sort!(deletion_candidates, rev = true)
        cand = map(x -> (x[1], kmerMatches[x[2]]), deletion_candidates)
        prev_remaining_matches = remaining_matches
        for j in 1 : max_pruning_num
            if deletion_candidates[j][1] > corScore && remaining_matches > 3
                kmer = deletion_candidates[j][2]
                deletionFlags[kmer] = true
                kmerMetrics = RemoveVector(kmerMetrics, kmerMatches[kmer].posA, kmerMatches[kmer].posB)
                corScore = CorrelationScore(kmerMetrics)
                remaining_matches -= 1
            else
                break
            end
        end
        if remaining_matches == prev_remaining_matches
            break
        end
    end

    #Ensure that kmers are compatible
    filteredKmerMatches = KmerMatch[]
    prevA = -k
    prevB = -k
    for i in 1 : matchCount
        match = kmerMatches[i]
        if !deletionFlags[i] && prevA + k <= match.posA && prevB + k <= match.posB
            prevA = match.posA
            prevB = match.posB
            push!(filteredKmerMatches, match)
        end
    end
    
    return filteredKmerMatches
end