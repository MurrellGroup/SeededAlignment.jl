# TODO test the function works
function levenshtein(seq1::BioSequence{T}, seq2::BioSequence{T}) where T<:Alphabet
    len1 = length(seq1)
    len2 = length(seq2)

    # Initialize distance matrix
    d = Array{Int}(undef, len1+1, len2+1)
    for i in 1:len1+1
        d[i, 1] = i-1
    end
    for j in 1:len2+1
        d[1, j] = j-1
    end

    # Compute distances
    for i in 2:len1+1
        for j in 2:len2+1
            cost = seq1[i-1] == seq2[j-1] ? 0 : 1
            d[i, j] = minimum((
                d[i-1, j] + 1,     # deletion
                d[i, j-1] + 1,     # insertion
                d[i-1, j-1] + cost # substitution
            ))
        end
    end

    return d[len1+1, len2+1]
end

function SP_score(target_alignment::Tuple{BioSequence{T},BioSequence{T}}, src_alignment::Tuple{BioSequence{T},BioSequence{T}}) where T<:Alphabet
    @assert length(target_alignment[1]) == length(src_alignment[1]) # FIXME handle case where different lengths if it arises enough
    seqlen = length(target_alignment[1])
    display_offset = 10
    matches = 0
    for j in 1:seqlen
        if target_alignment[1][j] == src_alignment[1][j] && target_alignment[2][j] == src_alignment[2][j]
            matches += 1
        else
            println("alignment_index: ", j)
            if j-display_offset > 0 && seqlen > j+display_offset
                println("target alignment:\n",target_alignment[1][j-display_offset:j+display_offset],"\n",target_alignment[2][j-display_offset:j+display_offset])
                println("src alignment:\n", src_alignment[1][j-display_offset:j+display_offset],"\n",src_alignment[2][j-display_offset:j+display_offset])
            elseif seqlen < j+display_offset
                println("target alignment:\n",target_alignment[1][j-display_offset:end],"\n",target_alignment[2][j-display_offset:end])
                println("src alignment:\n", src_alignment[1][j-display_offset:end],"\n",src_alignment[2][j-display_offset:end])
            else
                println("target alignment:\n",target_alignment[1][1:j+display_offset],"\n",target_alignment[2][1:j+display_offset])
                println("src alignment:\n", src_alignment[1][1:j+display_offset],"\n",src_alignment[2][1:j+display_offset])
            end
        end
    end
    return (matches/seqlen)
end

#= implement these on msa_level
function SP_score()
    seqlen = length()
end

function TC_score()

end

=#