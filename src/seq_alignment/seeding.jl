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
    min_match_score = maximum(t -> match_score_matrix[t, t], 1 : 4)
    gap_score_estimate = maximum(move -> move.score / move.step_length, (vgap_moves..., hgap_moves...))
    # if gap_score_estimate less than extension_score we take an average
    if gap_score_estimate < extension_score
        gap_score_estimate = (gap_score_estimate + 2 * extension_score) / 3
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

# helper types and functions
struct Endpoint
    x::Int64
    y::Int64
    id::Int64
    isBeginning::Bool
end

struct SampleMetrics{T<:Real}
    n::Int
    meanX::T
    meanY::T
    nVarX::T
    nVarY::T
    dotXY::T
    correlation::T
end

mutable struct KmerMatch
    posA::Int64
    posB::Int64
end

function CorrelationScore(s::SampleMetrics)
    return sqrt(s.n) * s.correlation
end

# Estimate pealty between kmer matches x and y
function getConnectionScore(x::KmerMatch, y::KmerMatch, k::Int64, gap_score_estimate::Float64, match_score_estimate::Float64)
    m = y.posA - x.posA - k
    n = y.posB - x.posB - k
    return abs(m - n) * gap_score_estimate + min(m, n) * match_score_estimate
end

function SampleMetrics(X::Vector{<:Integer}, Y::Vector{<:Integer})
    length(X) == length(Y) || error("Vectors must be of equal length.")
    n = length(X)
    meanX = sum(X) / n
    meanY = sum(Y) / n
    nVarX = sum(abs2, X .- meanX) #matchCount * variance of A
    nVarY = sum(abs2, Y .- meanY) #matchCount * variance of B
    dotXY = dot(X, Y)
    correlation = (dotXY - n * meanX * meanY) / sqrt(nVarX * nVarY)
    return SampleMetrics(n, meanX, meanY, nVarX, nVarY, dotXY, correlation)
end

function remove_point(s::SampleMetrics, x::Integer, y::Integer)
    n = s.n - 1
    meanX = (s.n * s.meanX - x) / n
    meanY = (s.n * s.meanY - y) / n
    nVarX = s.nVarX - (s.n / n) * (x - s.meanX)^2 #newMatchCount * new variance of A
    nVarY = s.nVarY - (s.n / n) * (y - s.meanY)^2 #newMatchCount * new variance of B
    dotXY = s.dotXY - x*y
    correlation = (dotXY - n * meanX * meanY) / sqrt(nVarX * nVarY)
    return SampleMetrics(n, meanX, meanY, nVarX, nVarY, dotXY, correlation)
end