using LinearAlgebra

mutable struct SampleMetrics2D
    n::Int64
    meanX::Float64
    meanY::Float64
    dotXY::Float64
    nVarX::Float64
    nVarY::Float64
    correlation::Float64
    SampleMetrics2D() = new()
end

function GetSampleMetrics2D(X::Array{Int}, Y::Array{Int})
    if length(X) != length(Y)
        error("Vectors must be of equal length.")
    end
    s = SampleMetrics2D()
    s.n = length(X)
    s.meanX = sum(X) / s.n
    s.meanY = sum(Y) / s.n
    s.dotXY = dot(X, Y)
    s.nVarX = sum(map(x -> (x - s.meanX)^2, X)) #matchCount * variance of A
    s.nVarY = sum(map(y -> (y - s.meanY)^2, Y)) #matchCount * variance of B
    s.correlation = (s.dotXY - s.n * s.meanX * s.meanY) / sqrt(s.nVarX * s.nVarY)
    return s
end

function RemoveVector(s::SampleMetrics2D, x, y)
    modS = SampleMetrics2D()
    modS.n = s.n - 1
    modS.meanX = (s.n * s.meanX - x) / modS.n
    modS.meanY = (s.n * s.meanY - y) / modS.n
    modS.dotXY = s.dotXY - x*y
    modS.nVarX = s.nVarX - (s.n / modS.n) * (x - s.meanX)^2 #newMatchCount * new variance of A
    modS.nVarY = s.nVarY - (s.n / modS.n) * (y - s.meanY)^2 #newMatchCount * new variance of B
    modS.correlation = (modS.dotXY - modS.n * modS.meanX * modS.meanY) / sqrt(modS.nVarX * modS.nVarY)
    return modS
end