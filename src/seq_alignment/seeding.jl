# helper types and functions
struct KmerMatch
    posA::Int64
    posB::Int64
end

# seeding functions

# Helper function and types for find_kmer_matches
# mutable struct optimized for kmerDict
mutable struct KmerPositions
    p1::UInt16
    p2::UInt16
    p3::UInt16
    p4::UInt16
    p5::UInt16
    skip::Bool
end

# Helper: insert a new position into the first free slot
@inline function insert!(kp::KmerPositions, pos::UInt16)
    if kp.p1 == UInt16(0)
        kp.p1 = pos
    elseif kp.p2 == UInt16(0)
        kp.p2 = pos
    elseif kp.p3 == UInt16(0)
        kp.p3 = pos
    elseif kp.p4 == UInt16(0)
        kp.p4 = pos
    elseif kp.p5 == UInt16(0)
        kp.p5 = pos
    else
        # purely to set that kmer should be ignored
        kp.skip = true
    end
end
const repetition_threshold::Int64 = 5

# finds all kmer_matches to be considered
@inbounds function find_kmer_matches(A::LongDNA{4}, B::LongDNA{4}, kmer_length::Int64; A_is_ref::Bool=false) 
    # Abbreviations
    k = kmer_length
    m = length(A)
    n = length(B)

    # List of all kmer matches between A and B
    kmer_matches = KmerMatch[]
    # Produce dictionary of all kmers in A
    kmerDict = Dict{UInt64, KmerPositions}()
    if !A_is_ref
        for i in 1:m-k+1
            hash_key = encode_kmer(A, i, k)
            if haskey(kmerDict, hash_key)
                KP = kmerDict[hash_key]
                if KP.skip == false
                    insert!(KP, i)
                end
            else
                kmerDict[hash_key] = KmerPositions(i,0,0,0,0, false)
            end
        end
    else
        # seeds in A must obide by reading frame
        for i in 1:(m-k+1)÷3
            hash_key = encode_kmer(A, 3*(i-1)+1, k)
            pos = UInt16(3*(i-1)+1)
            if haskey(kmerDict, hash_key)
                KP = kmerDict[hash_key]
                if KP.skip == false
                    insert!(KP, pos)
                end
            else
                kmerDict[hash_key] = KmerPositions(pos, UInt16(0), UInt16(0), UInt16(0), UInt16(0), false)
            end
        end
    end

    # Search B for any matching kmers
    diagonals = fill(typemin(UInt32), m+n)
    #= diagonals[i] = the rightmost kmer start_index in A that matches with a kmer of B in diagonal i. 
    Used to avoid overlapping kmerMatches=#

    #= NOTE: improvement idea if alignment quality is lacking - 
        IDEA: find all kmerMatches on diagional i. Sort them, itterate, check if current KmerMatch overlaps with previously added KmerMatch. 
        On the other hand this causes more kmerMatches to appear and make the alignment slower overall.
    =#
    for iB in 1 : n-k+1
        hash_key = encode_kmer(B, iB, k)
        # skips kmers which are too common in sequence A
        if (haskey(kmerDict, hash_key) && kmerDict[hash_key].skip == false)
            #Add found match(es) to list
            KP = kmerDict[hash_key]
            for iA in (KP.p1, KP.p2, KP.p3, KP.p4, KP.p5)
                iA == 0 && break
                diag_idx = iA - iB + n + 1
                #= 
                    NOTE: if iA is big early then it could impact alignment quality because it blocks 
                    other matches along that diagonal. Remedy was suggested above. But repetition_threshold makes this less likely
                =#
                # sufficent condition for no overlap
                if diagonals[diag_idx] + k <= iA
                    push!(kmer_matches, KmerMatch(iA, iB))
                    diagonals[diag_idx] = iA
                end
            end
        end
    end
    
    return kmer_matches
end

# structs and helper methods for kmer selection using approximated nw_align score
struct Endpoint
    x::Int64
    y::Int64
    id::Int64
    isBeginning::Bool
end
# bitmask for kmer hashing
function encode_kmer(A::LongDNA{4}, start::Int64, k::Int64)
    code::UInt64 = 0
    @inbounds @simd for j in 0:k-1
        code <<= 2
        code |= trailing_zeros(BioSequences.compatbits(A[start+j]))  # .data is 0,1,2,3 for A,C,G,T
    end
    return code
end
# Estimate pealty between kmer matches x and y
@inline function getConnectionScore(x::KmerMatch, y::KmerMatch, k::Int64, gap_score_estimate::Float64, match_score_estimate::Float64)
    m = y.posA - x.posA - k
    n = y.posB - x.posB - k
    return abs(m - n) * gap_score_estimate + min(m, n) * match_score_estimate
end

@inline function get_match_score_avg_and_min(M::Matrix{Float64}, match_score::Float64, mismatch_score::Float64)
    return sum(M)/16, maximum(t -> match_score_matrix[t, t], 1 : 4)
end

@inline function get_match_score_avg_and_min(::Nothing, match_score::Float64, mismatch_score::Float64)
    return match_score/4+(3/4)*mismatch_score, match_score
end
#kmer selection based on approximated nw_align score
function select_kmer_path(kmerMatches::Vector{KmerMatch}, m::Int64, n::Int64, match_score_matrix::S, match_score::Float64, mismatch_score::Float64,
    vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, extension_score::Float64, k::Int64) where {X, Y, S<:Union{Nothing, Matrix{Float64}}}
    
    # Produce constants used for estimating scores without A and B
    avg_match_score, min_match_score = get_match_score_avg_and_min(match_score_matrix, match_score, mismatch_score)
    match_score_estimate = avg_match_score - min_match_score
    gap_score_estimate = maximum(move -> move.score / move.step_length, (vgap_moves..., hgap_moves...))
    # if gap_score_estimate less than extension_score we take an average
    if gap_score_estimate < extension_score
        gap_score_estimate = (gap_score_estimate + 2 * extension_score) / 3
    end
    gap_score_estimate -= min_match_score/2
    matchCount = length(kmerMatches)
    
    #Produce list of endpoints - two for each kmer
    endpoints = Vector{Endpoint}(undef,matchCount * 2)
    for i in 1 : matchCount
        match = kmerMatches[i] 
        endpoints[2 * i - 1] = Endpoint(match.posA, match.posB, i, true)
        endpoints[2 * i] = Endpoint(match.posA + k, match.posB + k, i, false)
    end
    sort!(endpoints, by = e -> (e.x, e.y, e.isBeginning))
    
    # Transitions
    bestConnections = Vector{Int64}(undef, matchCount)
    bestScores = Vector{Float64}(undef, matchCount)

    #Chains speed up the DP
    chainNum = 0
    chains = Vector{Vector{Int64}}() # Disjoint sequences of compatible kmers
    chainIds = Vector{Int64}(undef, matchCount) # The chain of each kmer
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
                    if score > bestScores[e.id]
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
        kmer = argmax([bestScores[i] + getConnectionScore(kmerMatches[i], KmerMatch(m, n), k, gap_score_estimate, match_score_estimate) for i in 1 : matchCount])
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

# structs and helper methods for kmer selection based on correlation scoring
struct SampleMetrics
    n::Int
    meanX::Float64
    meanY::Float64
    nVarX::Float64
    nVarY::Float64
    dotXY::Float64 # this seems to always be Int64
    correlation::Float64
end

function SampleMetrics(X::Vector{Int64}, Y::Vector{Int64})
    length(X) == length(Y) || error("Vectors must be of equal length.")
    n = length(X)
    meanX = sum(X) / n
    meanY = sum(Y) / n
    nVarX = sum(abs2, X .- meanX) #matchCount * variance of A
    nVarY = sum(abs2, Y .- meanY) #matchCount * variance of B
    dotXY = 0
    @inbounds @simd for i in 1:length(X)
        dotXY += X[i] * Y[i]
    end
    correlation = (dotXY - n * meanX * meanY) / sqrt(nVarX * nVarY)
    return SampleMetrics(n, meanX, meanY, nVarX, nVarY, Float64(dotXY), correlation)
end

@inline function CorrelationScore(s::SampleMetrics)
    return sqrt(s.n) * s.correlation
end

# TODO benchmark if these benefit from inline
function remove_point(s::SampleMetrics, x::Int64, y::Int64)
    n = s.n - 1
    meanX = (s.n * s.meanX - x) / n
    meanY = (s.n * s.meanY - y) / n
    nVarX = s.nVarX - (s.n / n) * (x - s.meanX)^2 #newMatchCount * new variance of A
    nVarY = s.nVarY - (s.n / n) * (y - s.meanY)^2 #newMatchCount * new variance of B
    dotXY = s.dotXY - x*y
    correlation = (dotXY - n * meanX * meanY) / sqrt(nVarX * nVarY)
    return SampleMetrics(n, meanX, meanY, nVarX, nVarY, dotXY, correlation)
end

#Kmer selection variant algorithm based on correlation scoring
function select_max_correlation_kmer_path(kmerMatches::Vector{KmerMatch}, k::Int64)

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