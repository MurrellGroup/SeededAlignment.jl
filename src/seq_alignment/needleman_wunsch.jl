# Convert NucleicAcid to integer A -> 1, C -> 2, G -> 3, T -> 4
toInt(x::NucleicAcid) = trailing_zeros(reinterpret(UInt8, x)) + 1

function simple_match_penalty_matrix(match_score, mismatch_score, n=4)
    m = fill(mismatch_score, n, n)
    m[diagind(m)] .= match_score
    return m
end

const codon_table::Vector{AminoAcid} = AminoAcid[
    # AAA AAC AAG AAT
    AminoAcid('K'), AminoAcid('N'), AminoAcid('K'), AminoAcid('N'),
    # ACA ACC ACG ACT
    AminoAcid('T'), AminoAcid('T'), AminoAcid('T'), AminoAcid('T'),
    # AGA AGC AGG AGT
    AminoAcid('R'), AminoAcid('S'), AminoAcid('R'), AminoAcid('S'),
    # ATA ATC ATG ATT
    AminoAcid('I'), AminoAcid('I'), AminoAcid('M'), AminoAcid('I'),

    # CAA CAC CAG CAT
    AminoAcid('Q'), AminoAcid('H'), AminoAcid('Q'), AminoAcid('H'),
    # CCA CCC CCG CCT
    AminoAcid('P'), AminoAcid('P'), AminoAcid('P'), AminoAcid('P'),
    # CGA CGC CGG CGT
    AminoAcid('R'), AminoAcid('R'), AminoAcid('R'), AminoAcid('R'),
    # CTA CTC CTG CTT
    AminoAcid('L'), AminoAcid('L'), AminoAcid('L'), AminoAcid('L'),

    # GAA GAC GAG GAT
    AminoAcid('E'), AminoAcid('D'), AminoAcid('E'), AminoAcid('D'),
    # GCA GCC GCG GCT
    AminoAcid('A'), AminoAcid('A'), AminoAcid('A'), AminoAcid('A'),
    # GGA GGC GGG GGT
    AminoAcid('G'), AminoAcid('G'), AminoAcid('G'), AminoAcid('G'),
    # GTA GTC GTG GTT
    AminoAcid('V'), AminoAcid('V'), AminoAcid('V'), AminoAcid('V'),

    # TAA TAC TAG TAT
    AminoAcid('*'), AminoAcid('Y'), AminoAcid('*'), AminoAcid('Y'),
    # TCA TCC TCG TCT
    AminoAcid('S'), AminoAcid('S'), AminoAcid('S'), AminoAcid('S'),
    # TGA TGC TGG TGT
    AminoAcid('*'), AminoAcid('C'), AminoAcid('W'), AminoAcid('C'),
    # TTA TTC TTG TTT
    AminoAcid('L'), AminoAcid('F'), AminoAcid('L'), AminoAcid('F')
]

# TODO finish this fast implementation
function fast_translate(dna_seq::NTuple{3, DNA})
    # TODO unit test the codon_table
    hash_index = (toInt(dna_seq[1])-1)*16 + (toInt(dna_seq[2])-1)*4 + toInt(dna_seq[3])
    return codon_table[hash_index]
end

""" 
    nw_align(A::LongDNA{4},B::LongDNA{4},moveset::MoveSet,scoreScheme::ScoreScheme)
    
    Takes two ungapped LongDNA{4} sequences and computes an optimal pairwise alignment 
    with respect to the moveset and the scoreScheme. 

"""
# needleman_wunsch wrapper
function nw_align(A::LongDNA{4}, B::LongDNA{4}, m::MoveSet, s::ScoreScheme; clean_up_flag=false::Bool, codon_matching_enabled=false::Bool)
    # unpack arguments and call another nw_align
    nw_align(
        A, B, s.match_score, s.mismatch_score, m.match_moves,
        m.vert_moves, m.hor_moves, s.extension_score,
        s.edge_ext_begin, s.edge_ext_end, clean_up_flag,
        codon_matching_enabled, s.codon_match_bonus
    )
end

function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score::Float64, mismatch_score::Float64, 
                match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move})

    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), match_moves, vgap_moves, hgap_moves) 
end

#Needleman Wunsch alignment without affine scoring
function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2},
                match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move})

    n, m = length(A), length(B)

    #Margins in the dp_matrix eliminates branching by avoiding boundschecking
    column_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1
    row_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1

    #The matrix and the sequences are expanded according to the offsets
    column_boundary = n + column_offset
    row_boundary = m + row_offset
    
    #extend sequences to fix offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B

    #Initialize DP matrix
    #The cell at [x + column_offset, y + column_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64,row_boundary,column_boundary)
    # Assign score 0 to starting position
    dp_matrix[row_offset, column_offset] = 0.0
    # itterate through dp_matrix
    for row_index ∈ row_offset:row_boundary
        for column_index ∈ column_offset:column_boundary
            # calculate which positions of the sequences are being aligned. NOTE that the positions are 0 indexed
            top_sequence_pos = column_index-column_offset
            left_sequence_pos = row_index-row_offset
            # find the best diagonal move
            for k ∈ match_moves
                mismatch_sum = sum(t -> match_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : k.step)
                dp_matrix[row_index, column_index] = min(
                    dp_matrix[row_index, column_index],
                    dp_matrix[row_index-k.step,column_index-k.step]+k.score+mismatch_sum
                )
            end

            for k ∈ vgap_moves
                # the lowest score of vertical gaps is assigned to the current position in the matrix
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase && #TODO rename
                        (left_sequence_pos - k.step) % k.horizontal_stride == k.horizontal_phase && #TODO rename
                        dp_matrix[row_index-k.step, column_index] + k.score < dp_matrix[row_index, column_index]
                    # update current matrix position with lower score
                    dp_matrix[row_index, column_index] = dp_matrix[row_index - k.step, column_index] + k.score
                end
            end
            
            for k ∈ hgap_moves
                # the lowest score of horizontal gaps is asigned to the current position in the matrix if the score is lower than the current one
                if (top_sequence_pos - k.step) % k.vertical_stride == k.vertical_phase && #TODO rename
                        (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase && #TODO rename
                        dp_matrix[row_index, column_index - k.step] + k.score < dp_matrix[row_index, column_index]
                    # update current matrix position with lower score
                    dp_matrix[row_index, column_index] = dp_matrix[row_index, column_index - k.step] + k.score
                end
            end
        end
    end

    # Backtracking
    x = column_boundary
    y = row_boundary
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")
    while x > column_offset || y > row_offset
        top_sequence_pos = x-column_offset
        left_sequence_pos = y-row_offset
        if x == column_offset # first column
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == row_offset # first row
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # iterate through diagonal (match) moves
            for k ∈ match_moves
                # calculate total (mis-)match score from move
                s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                # check if the move leads to the current cell
                if dp_matrix[y,x] == dp_matrix[y - k.step,x - k.step] + k.score + s

                    # write the resulting sequences
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, B2[y - i])
                    end
                    x -= k.step
                    y -= k.step
                    break
                end
            end

            for k ∈ vgap_moves
                # iterate through possible vertical moves and check if the move leads to the current cell
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase && 
                        (left_sequence_pos - k.step) % k.horizontal_stride == k.horizontal_phase && 
                        dp_matrix[y, x] == dp_matrix[y- k.step,x] + k.score
                    # perform move
                    for i ∈ 1:k.step
                        push!(res_A, DNA_Gap)
                        push!(res_B, B2[y - i])
                    end
                    # update
                    y -= k.step
                    break
                end
            end

            for k ∈ hgap_moves
                # iterate through possible horizontal moves and check if move leads to the current cell
                if (top_sequence_pos - k.step) % k.vertical_stride == k.vertical_phase && 
                    (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase && 
                    dp_matrix[y,x] == dp_matrix[y,x-k.step] + k.score
                    #perform move
                    for i ∈ 1:k.step
                        push!(res_A, A2[x - i])
                        push!(res_B, DNA_Gap)
                    end
                    # update
                    x -= k.step
                    break
                end
            end

        end
    end

    return reverse(res_A), reverse(res_B)
end

# Needleman Wunsch alignment with affine scoring
function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score::Float64, mismatch_score::Float64, 
        match_moves::Vector{Move}, vgap_moves::Vector{Move}, hgap_moves::Vector{Move}, extension_score::Float64, 
        edge_extension_begin=false::Bool, edge_extension_end=false::Bool, clean_up_flag=false::Bool, 
        codon_matching_enabled=false::Bool, codon_match_score::Float64 = -1.0)
 
    nw_align(A, B, simple_match_penalty_matrix(match_score, mismatch_score), 
        match_moves, vgap_moves, hgap_moves, extension_score, edge_extension_begin, edge_extension_end, 
        clean_up_flag, codon_matching_enabled, codon_match_score) 
end

function nw_align(A::LongDNA{4}, B::LongDNA{4}, match_score_matrix::Array{Float64, 2}, match_moves::Vector{Move}, vgap_moves::Vector{Move},
    hgap_moves::Vector{Move}, extension_score::Float64, edge_extension_begin=false::Bool, edge_extension_end=false::Bool, 
    clean_up_flag=false::Bool, codon_matching_enabled=false::Bool, codon_match_bonus::Float64 =-1.0)
    
    n, m = length(A), length(B)

    # TODO get rid of this function call by including non-affine case in the affine one. 
    # TODO force or strongly recommend extension_score >= 0. 
    if !(extension_score > 0)
        # Do non-affine NW
        return nw_align(A, B, match_score_matrix, match_moves, vgap_moves, hgap_moves)
    end

    # Offset indicies to avoid bounds-checking
    column_offset = maximum(k -> k.step, vcat(match_moves, hgap_moves)) + 1
    row_offset = maximum(k -> k.step, vcat(match_moves, vgap_moves)) + 1
    column_boundary = n + column_offset
    row_boundary = m + row_offset

    # Length of sequences and matrices are increased according to offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B

    # Initialize DP matrix
    # The cell at [x + row_offset, y + column_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(Inf64, row_boundary, column_boundary)

    # Assign score 0 to the empty alignment
    dp_matrix[row_offset, column_offset] = 0.0

    # Affine moves requires two extra DP matrices
    vaffine_matrix = fill(Inf64, row_boundary, column_boundary)
    haffine_matrix = fill(Inf64, row_boundary, column_boundary)
    
    # allow starting in extending
    if edge_extension_begin
        vaffine_matrix[row_offset,column_offset] = 0.0
        haffine_matrix[row_offset,column_offset] = 0.0
    end

    # Main DP -step
    for row_index ∈ row_offset : row_boundary
        for column_index ∈ column_offset : column_boundary
            top_sequence_pos = column_index-column_offset
            left_sequence_pos = row_index-row_offset

            # reward for matching codons if enabled
            if codon_matching_enabled && (top_sequence_pos) % 3 == 0
                # performance optimized translation
                ref_AA = fast_translate((A2[column_index-3],A2[column_index-2],A2[column_index-1]))
                seq_AA = fast_translate((B2[row_index-3],B2[row_index-2],B2[row_index-1]))
                if  ref_AA == seq_AA
                    # reward codon_match
                    mismatch_sum = sum(t -> match_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : 3)
                    mismatch_sum += codon_match_bonus
                    # update dp_matrix
                    dp_matrix[row_index, column_index] = min(
                        dp_matrix[row_index, column_index], 
                        dp_matrix[row_index-3,column_index-3]+mismatch_sum
                    )
                end
            end

            # find the best diagonal move
            for k ∈ match_moves
                mismatch_sum = sum(t -> match_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : k.step)
                dp_matrix[row_index, column_index] = min(
                    dp_matrix[row_index, column_index], 
                    dp_matrix[row_index-k.step,column_index-k.step]+k.score+mismatch_sum
                )
            end

            # finds the best vertical move
            for k ∈ vgap_moves
                if (top_sequence_pos) % k.vertical_stride == k.vertical_phase &&
                    (left_sequence_pos-k.step) % k.horizontal_stride == k.horizontal_phase
                    if k.extensionAble # TODO rename
                        vaffine_matrix[row_index, column_index] = min(
                            vaffine_matrix[row_index, column_index],
                            vaffine_matrix[row_index - k.step, column_index] + extension_score * k.step,
                            dp_matrix[row_index - k.step, column_index] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index - k.step, column_index] + k.score
                        )
                    end
                end
            end

            # finds the best horizontal move
            for k ∈ hgap_moves
                if (top_sequence_pos-k.step) % k.vertical_stride == k.vertical_phase && 
                    (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase
                    if k.extensionAble # TODO rename
                        haffine_matrix[row_index, column_index] = min(
                            haffine_matrix[row_index, column_index],
                            haffine_matrix[row_index, column_index - k.step] + extension_score * k.step,
                            dp_matrix[row_index, column_index - k.step] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index, column_index - k.step] + k.score
                        )
                    end
                end
            end

            # find overall best move
            dp_matrix[row_index, column_index] = min(
                dp_matrix[row_index, column_index], 
                haffine_matrix[row_index, column_index],
                vaffine_matrix[row_index, column_index]
            )
        end
    end

    # handle ending alignment in extension state if enabled
    if edge_extension_end
        for row_index in row_offset : row_boundary
            vaffine_matrix[row_boundary,column_boundary] = min(
                vaffine_matrix[row_boundary,column_boundary],  
                dp_matrix[row_index,column_boundary] + extension_score*(row_boundary-row_index)
            )
        end
        for column_index in column_offset : column_boundary
            haffine_matrix[row_boundary,column_boundary] = min(
                haffine_matrix[row_boundary, column_boundary],
                dp_matrix[row_boundary,column_index] + extension_score*(column_boundary-column_index)
            )
        end
        # update end score
        dp_matrix[row_boundary, column_boundary] = min(
            dp_matrix[row_boundary, column_boundary], 
            haffine_matrix[row_boundary, column_boundary],
            vaffine_matrix[row_boundary, column_boundary]
        )
    end
   
    # Backtracking
    res_A = LongDNA{4}("")
    res_B = LongDNA{4}("")
    # Start at the final cell
    x = column_boundary
    y = row_boundary
    # Flags for affine moves
    must_move_ver = false
    must_move_hor = false
    while x > column_offset || y > row_offset
        top_sequence_pos = x-column_offset
        left_sequence_pos = y-row_offset
        if x == column_offset # first column
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == row_offset # first row
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # end extension backtrack
            if x == column_boundary && edge_extension_end
                for i in 1:y-row_offset
                    if isapprox(dp_matrix[y,x],dp_matrix[y-i,x]+i*extension_score)
                        for j ∈ 1 : i
                            push!(res_A, DNA_Gap)
                            push!(res_B, B2[y - j])
                        end
                        y -= i
                        # do stuff
                        break
                    end
                end
            end
            if y == row_boundary && edge_extension_end
                for i in 1:x-column_offset
                    if isapprox(dp_matrix[y,x],dp_matrix[y,x-i]+i*extension_score)
                        for j ∈ 1:i
                            push!(res_A, A2[x - j])
                            push!(res_B, DNA_Gap)
                        end
                        x -= i
                        # do stuff
                        break
                    end
                end
            end

            # record previous position
            px = x
            py = y

            # iterate through digonal match moves
            match_Found = false
            if !must_move_hor && !must_move_ver
                # reward for matching codons if enabled
                if codon_matching_enabled && (top_sequence_pos) % 3 == 0
                    mismatch_sum = sum(t -> match_score_matrix[toInt(A2[x - t]), toInt(B2[y - t])], 1 : 3)
                    ref_AA = fast_translate((A2[x-3],A2[x-2],A2[x-1]))
                    seq_AA = fast_translate((B2[y-3],B2[y-2],B2[y-1]))
                    if ref_AA == seq_AA
                        # reward codon_match
                        mismatch_sum += codon_match_bonus
                        # check if the move leads to the current cell
                        if isapprox(dp_matrix[y, x],dp_matrix[y - 3,x - 3] + mismatch_sum)
                            # record the path
                            for i ∈ 1:3
                                push!(res_A, A2[x - i])
                                push!(res_B, B2[y - i])
                            end
                            x -= 3
                            y -= 3
                            match_Found = true
                        end
                    end
                end

                if !match_Found
                    # go through matching moves
                    for k ∈ match_moves
                        # calculate total (mis-)match score
                        s = sum(t -> match_score_matrix[toInt(A2[x-t]), toInt(B2[y-t])], 1 : k.step)
                        # check if the move leads to the current cell
                        if isapprox(dp_matrix[y, x],dp_matrix[y - k.step,x - k.step] + k.score + s)
                            # record the path
                            for i ∈ 1:k.step
                                push!(res_A, A2[x - i])
                                push!(res_B, B2[y - i])
                            end
                            x -= k.step
                            y -= k.step
                            match_Found = true
                            break
                        end
                    end
                end
            end

            if !must_move_hor && !match_Found
                
                for k ∈ vgap_moves
                    !( ((top_sequence_pos) % k.vertical_stride == k.vertical_phase) &&
                       ((left_sequence_pos-k.step) % k.horizontal_stride == k.horizontal_phase) ) ? continue : 
                    # check if the move leads to the current cell
                    if k.extensionAble
                        current_score = must_move_ver ? vaffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (isapprox(current_score,vaffine_matrix[y-k.step, x] + extension_score * k.step))
                        can_move_regular = (isapprox(current_score,dp_matrix[y-k.step, x] + k.score))
                    else
                        current_score = dp_matrix[y,x]
                        can_move_affine = (false)
                        can_move_regular = (isapprox(current_score,dp_matrix[y-k.step,x] + k.score))
                    end
                    
                    if can_move_affine || can_move_regular
                        for i ∈ 1 : k.step
                            push!(res_A, DNA_Gap)
                            push!(res_B, B2[y - i])
                        end
                        y -= k.step

                        # constrain next move 
                        if !(y == row_offset)
                            must_move_ver = !can_move_regular
                        end
                        break
                    end
                end
            end

            if !must_move_ver && !match_Found

                # iterate through horizontal Move moves
                for k ∈ hgap_moves
                    !((top_sequence_pos-k.step) % k.vertical_stride == k.vertical_phase && 
                      (left_sequence_pos) % k.horizontal_stride == k.horizontal_phase) ? continue :
                    # check if the move leads to the current cell
                    if k.extensionAble
                        current_score = must_move_hor ? haffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (isapprox(current_score,haffine_matrix[y, x-k.step] + extension_score * k.step))
                        can_move_regular = (isapprox(current_score,dp_matrix[y,x-k.step] + k.score))
                    else
                        current_score = dp_matrix[y, x]
                        can_move_affine = (false)
                        can_move_regular = (isapprox(current_score,dp_matrix[y,x-k.step] + k.score))
                    end

                    if can_move_affine || can_move_regular
                        
                        for i ∈ 1:k.step
                            push!(res_A, A2[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step

                        # constrain next move
                        if !(x == column_offset)
                            must_move_hor = !can_move_regular
                        end
                        break
                    end
                end
            end

            # if no move was found
            if px == x && py == y
                error("Backtracking failed")
            end
        end
    end
    # full alignment 
    aligned_A = reverse(res_A)
    aligned_B = reverse(res_B)
    # clean_up single indels if enabled
    if clean_up_flag
        aligned_A, aligned_B = clean_alignment_readingframe(aligned_A, aligned_B)
    end
    # return alignment
    return aligned_A, aligned_B
end