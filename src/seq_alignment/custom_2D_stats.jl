struct SampleMetrics{T<:Real}
    n::Int
    meanX::T
    meanY::T
    nVarX::T
    nVarY::T
    dotXY::T
    correlation::T
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