""" 
	nw_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset=STD_NOISY_MOVESET, scoring::ScoringScheme=STD_SCORING)

(Needleman-wunsch wrapper - DE-NOVO)

Computes an optimal global pairwise alignment of the two ungapped DNA sequences `A` and `B`. This is done purely semantically without any awareness of protein encoding.

# Extended Help

# Arguments

- `A::LongDNA{4}`: 1st DNA sequence to be aligned
- `B::LongDNA{4}`: 2nd DNA sequence to be aligned
- `moveset::Moveset=STD_NOISY_MOVESET`: Defines allowable alignment moves (e.g. insertions/deletions and their penalty)
- `scoring::ScoreScheme=STD_SCORING`: Defines alignment scoring together with moveset

# Returns

- `Tuple{LongDNA{4},LongDNA{4}}`: Tuple representation of pairwise alignment of DNA sequences `A` and `B`.

# Examples
```julia
A = LongDNA{4}("AATGCTC")
B = LongDNA{4}("ACATGTC")
# produce alignment 
alignment = nw_align(A, B)
println(alignment)
#= resulting alignment:

alignment = (
	LongDNA{4}("A-ATGCTC"), 
	LongDNA{4}("ACATG-TC")
)

=#
```
"""
function nw_align(A::LongDNA{4}, B::LongDNA{4}; moveset::Moveset = STD_NOISY_MOVESET, scoring::ScoringScheme = STD_SCORING)   
    # throw exception if input sequences contains gaps
    !any(x -> x == DNA_Gap, A) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    !any(x -> x == DNA_Gap, B) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    # throw exception if invalid nucleotide letter in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), A) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), B) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
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
    _nw_align(
        A, B, moveset.vert_moves, moveset.hor_moves,
        scoring.extension_score, scoring.edge_ext_begin, scoring.edge_ext_end, scoring.nucleotide_match_score,
        scoring.nucleotide_mismatch_score, scoring.nucleotide_score_matrix, scoring.codon_match_bonus_score, codon_scoring_on, 
        do_clean_frameshifts, verbose
    )
end

""" 
	nw_align(; 
        	ref::LongDNA{4}, 
        	query::LongDNA{4}, 
        	moveset::Moveset = STD_CODON_MOVESET, 
        	scoring::ScoringScheme = STD_SCORING,
        	codon_scoring_on::Bool = true,
        	do_clean_frameshifts::Bool = false, 
        	verbose::Bool = false)

(Needleman-wunsch wrapper - CODING given trusted CDS anchor/reference)

Produces an optimal global pairwise alignment of two ungapped CDS (Coding DNA Sequences) `ref` and `query` by using an `ref` as an anchor to determine the apprioate reading frame coordinates for `query`. 

# Extended Help

# Arguments
- `ref::LongDNA{4}`: Anchored trusted CDS which decides the reading frame coordinates in the alignment
- `query::LongDNA{4}`: CDS (with possible frameshifts due to e.g. sequencing or annotation errors) which is aligned to `ref` and adopts its reading frame coordinates. 
- `moveset::Moveset = STD_CODON_MOVESET`: Defines allowable alignment moves (e.g. insertions/deletions and their penalty)
- `scoring::ScoringScheme = STD_SCORING`: Defines alignment scoring together with moveset
- `codon_scoring_on::Bool = true`: Whether to apply additional scoring on codon-level 
- `do_clean_frameshifts::Bool = false`: Whether to clean the alignment output of gaps which cause frameshifts (IMPORTANT: produces a protein alignment on a nucleotide level)
- `verbose::Bool = false`: Whether to verbosely display what edits were made during the cleaning of frameshifts. 

# Returns
- `Tuple{LongDNA{4},LongDNA{4}}`: Tuple representation of pairwise alignment of DNA sequences `ref` and `query`. 
Note that this represents a protein alignment on a nucleotide level if `(do_clean_frameshifts == true)`.

# Examples
```julia
anchor_CDS =    LongDNA{4}("ATGCCAGTA")
# untrusted_CDS may contain some frameshift errors due to e.g. sequencing or annotation errors.
untrusted_CDS = LongDNA{4}("ATGTA") 
# frameshift errors are removed from the cleaned alignment
cleaned_CDS_alignment = nw_align(ref=anchor_CDS, query=untrusted_CDS, clean_frameshifts=true)
println(cleaned_CDS_alignment)
#= resulting alignment:

cleaned_CDS_alignment = (
	LongDNA{4}("ATGCCAGTA"), 
	LongDNA{4}("ATG---NTA")
)
Here 'N' denotes ambigious nucleotide.
=#
```
"""
function nw_align(; 
    ref::LongDNA{4}, 
    query::LongDNA{4}, 
    moveset::Moveset = STD_CODON_MOVESET, 
    scoring::ScoringScheme = STD_SCORING,
    codon_scoring_on::Bool = true,
    do_clean_frameshifts::Bool = false, 
    verbose::Bool = false)
    # throw exception if input sequences contains gaps
    !any(x -> x == DNA_Gap, ref) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    !any(x -> x == DNA_Gap, query) || throw(ArgumentError("Input sequence contains gaps! - Sequences must be ungapped!"))
    # throw exception if invalid nucleotide letter in LongDNA{4}
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), ref) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), query) || throw(ArgumentError("Input sequence contains non-standard nucleotides! \nThe only accepted symbols are 'A', 'C', 'T' and 'G'"))
    # check that moveset takes reference reading frame into account
    contains_ref_move(moveset) || throw(ArgumentError("Invalid Moveset for reference to query alignment!\n", 
                                        "At least one Move in Moveset must consider reference reading (Move.ref=true)",
                                        " - in other words codon insertions or deletions must be allowed."))
    # unpack arguments and call the internal alignment function
    _nw_align(
        ref, query, moveset.vert_moves, moveset.hor_moves,
        scoring.extension_score, scoring.edge_ext_begin, scoring.edge_ext_end, scoring.nucleotide_match_score,
        scoring.nucleotide_mismatch_score, scoring.nucleotide_score_matrix, scoring.codon_match_bonus_score, codon_scoring_on, 
        do_clean_frameshifts, verbose
    )
end

# Needleman Wunsch alignment with affine scoring (internal function)
@inbounds @fastmath function _nw_align(A::LongDNA{4}, B::LongDNA{4}, vgap_moves::NTuple{X,Move}, hgap_moves::NTuple{Y,Move}, 
    extension_score::Float64=-0.3, edge_extension_begin=false::Bool, edge_extension_end=false::Bool, nuc_match_score::Float64=0.0,
    nuc_mismatch_score::Float64=-0.8, nuc_score_matrix::S=nothing, codon_match_bonus_score::Float64=6.0, codon_scoring_on=false::Bool, 
    do_clean_frameshifts=false::Bool, verbose=false::Bool) where {X, Y, S<:Union{Nothing, Matrix{Float64}}}

    n, m = length(A), length(B)

    # Offset indicies to avoid bounds-checking
    column_offset = maximum(k -> k.step_length, hgap_moves) + 1
    row_offset =    maximum(k -> k.step_length, vgap_moves) + 1
    column_boundary = n + column_offset
    row_boundary = m + row_offset

    # Length of sequences and matrices are increased according to offset
    A2 = LongDNA{4}("A")^(column_offset - 1) * A
    B2 = LongDNA{4}("A")^(row_offset - 1) * B
    # initialize stop_aa
    stop_aa = AminoAcid('*')
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
            # reward for matching codons if enabled
            if codon_scoring_on && (top_sequence_pos) % 3 == 0
                # performance optimized translation
                ref_AA = fast_translate(A2[column_index-3],A2[column_index-2],A2[column_index-1])
                seq_AA = fast_translate(B2[row_index-3],B2[row_index-2],B2[row_index-1])
                # check if AminoAcids match and handle stop_codons
                if ref_AA == seq_AA && ref_AA != stop_aa
                    # get nucleotide mismatch_score
                    match_score = sum(t -> score_match(A2[column_index-t], B2[row_index-t], nuc_score_matrix, nuc_match_score, nuc_mismatch_score),1:3)
                    # add codon match bonus
                    match_score += codon_match_bonus_score
                    # update dp_matrix
                    dp_matrix[row_index, column_index] = max(
                        dp_matrix[row_index, column_index],
                        dp_matrix[row_index-3,column_index-3]+match_score
                    ) 
                end
            end

            # check score if match current pair of nucleotides
            match_score = score_match(A2[column_index-1], B2[row_index-1], nuc_score_matrix, nuc_match_score, nuc_mismatch_score)
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
    # end extension backtrack
    if edge_extension_end && x == column_boundary 
        for i in 1:y-row_offset
            if fast_simpler_isapprox(dp_matrix[y,x],dp_matrix[y-i,x]+i*extension_score)
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
    # end extension backtrack
    if edge_extension_end && y == row_boundary
        for i in 1:x-column_offset
            if fast_simpler_isapprox(dp_matrix[y,x],dp_matrix[y,x-i]+i*extension_score)
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
    # loop rest of backtrack
    while x > column_offset || y > row_offset
        top_sequence_pos = x-column_offset
        if x == column_offset # first column
            push!(res_A, DNA_Gap)
            push!(res_B, B2[y - 1])
            y -= 1
        elseif y == row_offset # first row
            push!(res_A, A2[x - 1])
            push!(res_B, DNA_Gap)
            x -= 1
        else
            # record previous position
            px = x
            py = y
            # iterate through digonal match moves
            if !must_move_hor && !must_move_ver 
                # reward for matching codons if enabled
                if codon_scoring_on && (top_sequence_pos) % 3 == 0
                    ref_AA = fast_translate(A2[x-3],A2[x-2],A2[x-1])
                    seq_AA = fast_translate(B2[y-3],B2[y-2],B2[y-1])
                    if ref_AA == seq_AA && ref_AA != stop_aa
                        # get nucleotide mismatch_score
                        match_score = sum(t -> score_match(A2[x - t], B2[y - t], nuc_score_matrix, nuc_match_score, nuc_mismatch_score),1:3)
                        # add codon match bouns
                        match_score += codon_match_bonus_score
                        # check if the move leads to the current cell
                        if fast_simpler_isapprox(dp_matrix[y, x], dp_matrix[y - 3, x - 3] + match_score)
                            # record the path
                            for i ∈ 1:3
                                push!(res_A, A2[x - i])
                                push!(res_B, B2[y - i])
                            end
                            x -= 3
                            y -= 3
                            continue
                        end
                    end
                end

                # calculate total (mis-)match score
                s = score_match(A2[x-1], B2[y-1], nuc_score_matrix, nuc_match_score, nuc_mismatch_score)
                # check if the move leads to the current cell
                if fast_simpler_isapprox(dp_matrix[y, x],dp_matrix[y - 1, x - 1] + s)
                    # record the path
                    push!(res_A, A2[x - 1])
                    push!(res_B, B2[y - 1])
                    x -= 1
                    y -= 1
                    continue
                end
            end

            if !must_move_hor
                
                for k ∈ vgap_moves
                    !( !k.ref || (top_sequence_pos) % 3 == 0 ) ? continue : 
                    # check if the move leads to the current cell
                    if k.extendable
                        current_score = must_move_ver ? vaffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (fast_simpler_isapprox(current_score,vaffine_matrix[y-k.step_length, x] + extension_score * k.step_length))
                        can_move_regular = (fast_simpler_isapprox(current_score,dp_matrix[y-k.step_length, x] + k.score))
                    else
                        current_score = dp_matrix[y,x]
                        can_move_affine = (false)
                        can_move_regular = (fast_simpler_isapprox(current_score,dp_matrix[y-k.step_length,x] + k.score))
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

            if !must_move_ver

                # iterate through horizontal Move moves
                for k ∈ hgap_moves
                    !(!k.ref || (top_sequence_pos-k.step_length) % 3 == 0) ? continue :
                    # check if the move leads to the current cell
                    if k.extendable
                        current_score = must_move_hor ? haffine_matrix[y, x] : dp_matrix[y, x]
                        can_move_affine = (fast_simpler_isapprox(current_score,haffine_matrix[y, x-k.step_length] + extension_score * k.step_length))
                        can_move_regular = (fast_simpler_isapprox(current_score,dp_matrix[y,x-k.step_length] + k.score))
                    else
                        current_score = dp_matrix[y, x]
                        can_move_affine = (false)
                        can_move_regular = (fast_simpler_isapprox(current_score,dp_matrix[y,x-k.step_length] + k.score))
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
            if (px == x && py == y)
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
# specialized scoring functions
@inline score_match(a::DNA, b::DNA, ::Nothing, match_score::Float64, mismatch_score::Float64) = (a == b) ? match_score : mismatch_score
@inline score_match(a::DNA, b::DNA, M::Matrix{Float64}, match_score::Float64, mismatch_score::Float64) = M[toInt(a), toInt(b)]
# Didn't give noticiable performance boost
@inline function fast_simpler_isapprox(a::Float64, b::Float64; eps::Float64=1e-5)
    return abs(a - b) < eps
end