""" 
    
    nw_align(A::LongDNA{4},B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING
    
Needleman_Wunsch wrapper - no reference, i.e. makes no assumptions about the two sequences. 

Computes an optimal global pairwise alignment of two ungapped `LongDNA{4}` sequence  
with respect to the `Moveset` and the `ScoringScheme`. 

"""
function nw_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING)
    # check no reference informed moves (Moev.ref=true) in moveset
    !contains_ref_move(moveset) || throw(ArgumentError("Moveset contains move(s) that considers reading frame (Move.ref=true)", 
                                                       " when no reference sequence was given!\n","Either set Move.ref=false for all moves in moveset", 
                                                    "if you want to align without a reference sequence. Specify by nw_align(ref=seq1, query=seq2)."))
    # force no clean_up
    do_clean_frameshifts=false
    verbose=false
    # force no codon_scoring_on
    codon_scoring_on=false
    # unpack arguments and call the internal alignment function
    nw_align(
        A, B, moveset.vert_moves, moveset.hor_moves, 
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_score_matrix,
        scoring.edge_ext_begin, scoring.edge_ext_end, codon_scoring_on, do_clean_frameshifts, verbose
    )
end

""" 
    
    nw_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING,
        do_clean_frameshifts=false::Bool, verbose=false::Bool, codon_scoring_on=true::Bool)

Needleman_Wunsch wrapper - reference informed, i.e. assumes one of the sequence has intact reading frame. 

Optimally aligns a `query` sequence to a `ref` sequence using a codon-aware `Moveset` and `ScoringScheme`.

**NOTE** We always assume the readingFrame is 1
"""
function nw_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING, 
    do_clean_frameshifts=false::Bool, verbose=false::Bool, codon_scoring_on=true::Bool)

    # check that moveset takes reference reading frame into account
    contains_ref_move(moveset) || throw(ArgumentError("Invalid Moveset for reference to query alignment!\n", 
                                        "At least one Move in Moveset must consider reference reading (Move.ref=true)",
                                        " - in other words codon insertions or deletions must be allowed."))
    # unpack arguments and call the internal alignment function
    nw_align(
        ref, query, moveset.vert_moves, moveset.hor_moves, 
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_score_matrix,
        scoring.edge_ext_begin, scoring.edge_ext_end, codon_scoring_on, do_clean_frameshifts, verbose
    )
end

# Needleman Wunsch alignment with affine scoring (internal function)
function nw_align(A::LongDNA{4}, B::LongDNA{4}, vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, 
    nuc_score_matrix::Matrix{Float64}, extension_score::Float64, codon_score_matrix::Matrix{Float64}=BLOSUM62, edge_extension_begin=false::Bool, 
    edge_extension_end=false::Bool, codon_scoring_on=false::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool) where {X, Y}
    
    # throw exception if invalid alphabet in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))

    n, m = length(A), length(B)

    # Offset indicies to avoid bounds-checking
    column_offset = maximum(k -> k.step_length, hgap_moves) + 1
    row_offset =    maximum(k -> k.step_length, vgap_moves) + 1
    column_boundary = n + column_offset
    row_boundary = m + row_offset

    # Length of sequences and matrices are increased according to offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B

    # Initialize DP matrix
    # The cell at [x + row_offset, y + column_offset] is the score of the best alignment of A[1 : x] with B[1 : y]
    dp_matrix = fill(-Inf64, row_boundary, column_boundary)

    # Assign score 0 to the empty alignment
    dp_matrix[row_offset, column_offset] = 0.0

    # Affine moves requires two extra DP matrices
    vaffine_matrix = fill(-Inf64, row_boundary, column_boundary)
    haffine_matrix = fill(-Inf64, row_boundary, column_boundary)
    
    # allow starting in extending
    if edge_extension_begin
        vaffine_matrix[row_offset,column_offset] = 0.0
        haffine_matrix[row_offset,column_offset] = 0.0
    end

    # Main DP -step
    for row_index ∈ row_offset : row_boundary
        for column_index ∈ column_offset : column_boundary
            top_sequence_pos = (column_index-column_offset)
            left_sequence_pos = (row_index-row_offset)
            # reward for matching codons if enabled
            if codon_scoring_on && (top_sequence_pos) % 3 == 0
                # performance optimized translation
                ref_AA = fast_translate((A2[column_index-3],A2[column_index-2],A2[column_index-1]))
                seq_AA = fast_translate((B2[row_index-3],B2[row_index-2],B2[row_index-1]))
                # TODO handle stop codon
                if ref_AA != AminoAcid('*') && seq_AA != AminoAcid('*')
                    # get codon_mismatch_score
                    codon_match_score = codon_score_matrix[toInt(ref_AA),toInt(seq_AA)]
                    # get nucleotide mismatch_score
                    match_score = sum(t -> nuc_score_matrix[toInt(A2[column_index - t]), toInt(B2[row_index - t])], 1 : 3)
                    match_score += codon_match_score
                    # update dp_matrix
                    dp_matrix[row_index, column_index] = max(
                        dp_matrix[row_index, column_index],
                        dp_matrix[row_index-3,column_index-3]+match_score
                    ) 
                end
            end

            # check score if match current pair of nucleotides
            match_score = nuc_score_matrix[toInt(A2[column_index - 1]), toInt(B2[row_index - 1])]
            dp_matrix[row_index, column_index] = max(
                dp_matrix[row_index, column_index], 
                dp_matrix[row_index-1,column_index-1]+match_score
            )

            # finds the best vertical move
            for k ∈ vgap_moves
                if !k.ref || (top_sequence_pos) % 3 == 0
                    if k.extendable
                        vaffine_matrix[row_index, column_index] = max(
                            vaffine_matrix[row_index, column_index],
                            vaffine_matrix[row_index - k.step_length, column_index] + extension_score * k.step_length,
                            dp_matrix[row_index - k.step_length, column_index] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = max(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index - k.step_length, column_index] + k.score
                        )
                    end
                end
            end

            # finds the best horizontal move
            for k ∈ hgap_moves
                if !k.ref || (top_sequence_pos-k.step_length) % 3 == 0
                    if k.extendable
                        haffine_matrix[row_index, column_index] = max(
                            haffine_matrix[row_index, column_index],
                            haffine_matrix[row_index, column_index - k.step_length] + extension_score * k.step_length,
                            dp_matrix[row_index, column_index - k.step_length] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = max(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index, column_index - k.step_length] + k.score
                        )
                    end
                end
            end

            # find overall best move
            dp_matrix[row_index, column_index] = max(
                dp_matrix[row_index, column_index], 
                haffine_matrix[row_index, column_index],
                vaffine_matrix[row_index, column_index]
            )
        end
    end

    # handle ending alignment in extension state if enabled
    if edge_extension_end
        for row_index in row_offset : row_boundary
            vaffine_matrix[row_boundary,column_boundary] = max(
                vaffine_matrix[row_boundary,column_boundary],  
                dp_matrix[row_index,column_boundary] + extension_score*(row_boundary-row_index)
            )
        end
        for column_index in column_offset : column_boundary
            haffine_matrix[row_boundary,column_boundary] = max(
                haffine_matrix[row_boundary, column_boundary],
                dp_matrix[row_boundary,column_index] + extension_score*(column_boundary-column_index)
            )
        end
        # update end score
        dp_matrix[row_boundary, column_boundary] = max(
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
            if !must_move_hor && !must_move_ver # TODO check if we ever get out of index range here
                # reward for matching codons if enabled
                if codon_scoring_on && (top_sequence_pos) % 3 == 0
                    ref_AA = fast_translate((A2[x-3],A2[x-2],A2[x-1]))
                    seq_AA = fast_translate((B2[y-3],B2[y-2],B2[y-1]))
                    # TODO handle stop codon more flexibly
                    if ref_AA != AminoAcid('*') && seq_AA != AminoAcid('*')
                        # get codon_match_score
                        match_score = codon_score_matrix[toInt(ref_AA), toInt(seq_AA)]
                        # get nucleotide_match_score
                        match_score += sum(t -> nuc_score_matrix[toInt(A2[x - t]), toInt(B2[y - t])], 1 : 3)
                        # check if the move leads to the current cell
                        if isapprox(dp_matrix[y, x],dp_matrix[y - 3,x - 3] + match_score)
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
                    # calculate total (mis-)match score
                    s = nuc_score_matrix[toInt(A2[x-1]), toInt(B2[y-1])]
                    # check if the move leads to the current cell
                    if isapprox(dp_matrix[y, x],dp_matrix[y - 1,x - 1] + s)
                        # record the path
                        push!(res_A, A2[x - 1])
                        push!(res_B, B2[y - 1])
                        x -= 1
                        y -= 1
                        match_Found = true
                    end
                end
            end

            if !must_move_hor && !match_Found
                
                for k ∈ vgap_moves
                    !( !k.ref || (top_sequence_pos) % 3 == 0 ) ? continue : 
                    # check if the move leads to the current cell
                    if k.extendable
                        current_score = must_move_ver ? vaffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (isapprox(current_score,vaffine_matrix[y-k.step_length, x] + extension_score * k.step_length))
                        can_move_regular = (isapprox(current_score,dp_matrix[y-k.step_length, x] + k.score))
                    else
                        current_score = dp_matrix[y,x]
                        can_move_affine = (false)
                        can_move_regular = (isapprox(current_score,dp_matrix[y-k.step_length,x] + k.score))
                    end
                    
                    if can_move_affine || can_move_regular
                        for i ∈ 1 : k.step_length
                            push!(res_A, DNA_Gap)
                            push!(res_B, B2[y - i])
                        end
                        y -= k.step_length

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
                    !(!k.ref || (top_sequence_pos-k.step_length) % 3 == 0) ? continue :
                    # check if the move leads to the current cell
                    if k.extendable
                        current_score = must_move_hor ? haffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (isapprox(current_score,haffine_matrix[y, x-k.step_length] + extension_score * k.step_length))
                        can_move_regular = (isapprox(current_score,dp_matrix[y,x-k.step_length] + k.score))
                    else
                        current_score = dp_matrix[y, x]
                        can_move_affine = (false)
                        can_move_regular = (isapprox(current_score,dp_matrix[y,x-k.step_length] + k.score))
                    end

                    if can_move_affine || can_move_regular
                        
                        for i ∈ 1:k.step_length
                            push!(res_A, A2[x - i])
                            push!(res_B, DNA_Gap)
                        end
                        x -= k.step_length

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
    if do_clean_frameshifts
        aligned_A, aligned_B = clean_frameshifts(aligned_A, aligned_B, verbose=verbose)
    end
    # return alignment
    return aligned_A, aligned_B
end