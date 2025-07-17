""" 
    
    nw_align(A::LongDNA{4},B::LongDNA{4}; moveset::Moveset=std_noisy_moveset(), scoring::ScoringScheme=std_scoring()
    
Needleman_Wunsch wrapper - no reference, i.e. makes no assumptions about the two sequences. 

Computes an optimal global pairwise alignment of two ungapped `LongDNA{4}` sequence  
with respect to the `Moveset` and the `ScoringScheme`. 

"""
function nw_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=std_noisy_moveset(), scoring::ScoringScheme=std_scoring())
    !contains_ref_move(moveset) || throw(ArgumentError("moveset contains move(s) that consider "))
    # force no clean_up
    do_clean_frameshifts=false
    verbose=false
    # force no match_codons
    match_codons=false
    # unpack arguments and call the internal alignment function
    nw_align(
        A, B, moveset.vert_moves, moveset.hor_moves, 
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_match_bonus,
        scoring.edge_ext_begin, scoring.edge_ext_end, match_codons, do_clean_frameshifts, verbose
    )
end

""" 
    
    nw_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=std_codon_moveset(), scoring::ScoringScheme=std_scoring(),
        do_clean_frameshifts=false::Bool, verbose=false::Bool, match_codons=true::Bool)

Needleman_Wunsch wrapper - reference informed, i.e. assumes one of the sequence has intact reading frame. 

Optimally aligns a `query` sequence to a `ref` sequence using a codon-aware `Moveset` and `ScoringScheme`.

**NOTE** We always assume the readingFrame is 1
"""
function nw_align(; ref::LongDNA{4}, query::LongDNA{4}, moveset::Moveset=std_codon_moveset(), scoring::ScoringScheme=std_scoring(), 
    do_clean_frameshifts=false::Bool, verbose=false::Bool, match_codons=true::Bool)

    # check that moveset takes reference reading frame into account
    contains_ref_move(moveset) || throw(ArgumentError("Invalid Moveset for reference to query alignment!\n", 
                                        "At least one Move in Moveset must consider reference reading (Move.ref=true)",
                                        " - in other words codon insertions or deletions must be allowed."))

    # unpack arguments and call the internal alignment function
    nw_align(
        ref, query, moveset.vert_moves, moveset.hor_moves, 
        scoring.nucleotide_score_matrix, scoring.extension_score, scoring.codon_match_bonus,
        scoring.edge_ext_begin, scoring.edge_ext_end, match_codons, do_clean_frameshifts, verbose
    )
end

# Needleman Wunsch alignment with affine scoring (internal function)
function nw_align(A::LongDNA{4}, B::LongDNA{4}, vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, 
    match_score_matrix::Matrix{Float64}, extension_score::Float64, codon_match_bonus::Float64 =-2.0, edge_extension_begin=false::Bool, 
    edge_extension_end=false::Bool, match_codons=false::Bool, do_clean_frameshifts=false::Bool, verbose=false::Bool) where {X, Y}
    
    # throw exception if invalid alphabet in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))

    n, m = length(A), length(B)
    # Do non-affine NW
    if !(extension_score > 0)
        println("non-affine mode")
        extension_score = Inf64
        edge_ext_begin = false
        edge_ext_end = false
    end

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
            top_sequence_pos = (column_index-column_offset)
            left_sequence_pos = (row_index-row_offset)
            # reward for matching codons if enabled
            if match_codons && (top_sequence_pos) % 3 == 0 # FIXME bug keeps occuring here for some reason. 
                # performance optimized translation
                ref_AA = fast_translate((A2[column_index-3],A2[column_index-2],A2[column_index-1])) # FIXME bug keeps occuring here for some reason.
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

            # check score if match current pair of nucleotides
            match_score = match_score_matrix[toInt(A2[column_index - 1]), toInt(B2[row_index - 1])]
            dp_matrix[row_index, column_index] = min(
                dp_matrix[row_index, column_index], 
                dp_matrix[row_index-1,column_index-1]+match_score
            )

            # finds the best vertical move
            for k ∈ vgap_moves
                if !k.ref || (top_sequence_pos) % 3 == 0
                    if k.extendable
                        vaffine_matrix[row_index, column_index] = min(
                            vaffine_matrix[row_index, column_index],
                            vaffine_matrix[row_index - k.step_length, column_index] + extension_score * k.step_length,
                            dp_matrix[row_index - k.step_length, column_index] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
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
                        haffine_matrix[row_index, column_index] = min(
                            haffine_matrix[row_index, column_index],
                            haffine_matrix[row_index, column_index - k.step_length] + extension_score * k.step_length,
                            dp_matrix[row_index, column_index - k.step_length] + k.score
                        )
                    else
                        dp_matrix[row_index,column_index] = min(
                            dp_matrix[row_index,column_index],
                            dp_matrix[row_index, column_index - k.step_length] + k.score
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
                if match_codons && (top_sequence_pos) % 3 == 0
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
                    # calculate total (mis-)match score
                    s = match_score_matrix[toInt(A2[x-1]), toInt(B2[y-1])]
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

# helper functions and constants

# Convert NucleicAcid to integer A -> 1, C -> 2, G -> 3, T -> 4
toInt(x::NucleicAcid) = trailing_zeros(reinterpret(UInt8, x)) + 1

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

function fast_translate(dna_seq::NTuple{3, DNA})
    # TODO unit test the codon_table
    hash_index = (toInt(dna_seq[1])-1)*16 + (toInt(dna_seq[2])-1)*4 + toInt(dna_seq[3])
    return codon_table[hash_index]
end