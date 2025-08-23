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

function SP_score(target_alignment::Tuple{BioSequence{T},BioSequence{T}}, src_alignment::Tuple{BioSequence{T},BioSequence{T}}; verbose::Bool=false) where T<:Alphabet
    @assert length(target_alignment[1]) == length(src_alignment[1])
    seqlen = length(target_alignment[1])
    display_offset = 20
    matches = 0
    last_err_index = -display_offset
    for j in 1:seqlen
        if target_alignment[1][j] == src_alignment[1][j] && target_alignment[2][j] == src_alignment[2][j]
            matches += 1
        elseif j-last_err_index > display_offset && verbose
            last_err_index = j
            println("comparision of source and target alignment - alignment_index: ", j)
            if j-display_offset > 0 && seqlen > j+display_offset
                println("target alignment:")
                color_diff(target_alignment[1][j-display_offset:j+display_offset],target_alignment[2][j-display_offset:j+display_offset])
                println("src alignment:") 
                color_diff(src_alignment[1][j-display_offset:j+display_offset],src_alignment[2][j-display_offset:j+display_offset])
                println("\nDifference between the 1st sequences")
                print("target alignment[1]: ")
                color_diff(target_alignment[1][j-display_offset:j+display_offset], src_alignment[1][j-display_offset:j+display_offset], second=false)
                print("src alignment[1]:    ") 
                color_diff(target_alignment[1][j-display_offset:j+display_offset], src_alignment[1][j-display_offset:j+display_offset], first=false)
                println("\nDifference between the 2nd sequences")
                print("target alignment[2]: ")
                color_diff(target_alignment[2][j-display_offset:j+display_offset], src_alignment[2][j-display_offset:j+display_offset], second=false)
                print("src alignment[2]:    ") 
                color_diff(target_alignment[2][j-display_offset:j+display_offset], src_alignment[2][j-display_offset:j+display_offset], first=false)
                println("--------------------------------------------------------------------------------------------------------------------------")
            elseif seqlen < j+display_offset
                println("target alignment:"),
                color_diff(target_alignment[1][j-display_offset:end],target_alignment[2][j-display_offset:end])
                println("src alignment:") 
                color_diff(src_alignment[1][j-display_offset:end],src_alignment[2][j-display_offset:end])
                println("\nDifference between 1st sequences")
                print("target alignment[1]: ")
                color_diff(target_alignment[1][j-display_offset:end], src_alignment[1][j-display_offset:end], second=false)
                print("src alignment[1]:    ") 
                color_diff(target_alignment[1][j-display_offset:end], src_alignment[1][j-display_offset:end], first=false)
                println("\nDifference between 2nd sequences")
                print("target alignment[2]: ")
                color_diff(target_alignment[2][j-display_offset:end], src_alignment[2][j-display_offset:end], second=false)
                print("src alignment[2]:    ") 
                color_diff(target_alignment[2][j-display_offset:end], src_alignment[2][j-display_offset:end], first=false)
                println("--------------------------------------------------------------------------------------------------------------------------")
            else
                println("target alignment:")
                color_diff(target_alignment[1][1:j+display_offset],target_alignment[2][1:j+display_offset])
                println("src alignment:") 
                color_diff(src_alignment[1][1:j+display_offset],src_alignment[2][1:j+display_offset])
                println("\nDifference between 1st sequences")
                print("target alignment[1]: ")
                color_diff(target_alignment[1][1:j+display_offset], src_alignment[1][1:j+display_offset], second=false)
                print("src alignment[1]:    ") 
                color_diff(target_alignment[1][1:j+display_offset], src_alignment[1][1:j+display_offset], first=false)
                println("\nDifference between 2nd sequences")
                print("target alignment[2]: ")
                color_diff(target_alignment[2][1:j+display_offset], src_alignment[2][1:j+display_offset], second=false)
                print("src alignment[2]:    ") 
                color_diff(target_alignment[2][1:j+display_offset], src_alignment[2][1:j+display_offset], first=false)
                println("--------------------------------------------------------------------------------------------------------------------------")
            end
        end
    end
    return (matches/seqlen)
end

function color_diff(s1::BioSequence{T}, s2::BioSequence{T}; first=true, second=true) where T<:Alphabet
    red   = "\033[31m"
    reset = "\033[0m"

    n = max(length(s1), length(s2))
    s1 = lpad(s1, n)
    s2 = lpad(s2, n)
    if first
        # First string
        for i in 1:n
            if s1[i] == s2[i]
                print(s1[i])               # normal
            else
                print(red, s1[i], reset)   # colored
            end
        end
        println()
    end

    if second
        # Second string
        for i in 1:n
            if s1[i] == s2[i]
                print(s2[i])               # normal
            else
                print(red, s2[i], reset)   # colored
            end
        end
        println()
    end
end