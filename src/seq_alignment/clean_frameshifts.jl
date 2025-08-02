
"""
    clean_frameshifts(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4}; verbose::Bool=false)

Takes a pairwise alignment of a reference (with known reading frame) and a sequence, and removes frameshift mutations 
which don't respect the reference's reading frame. This is done by removing insertions from the alignment or inserting 
ambigious nucleotides into deletions.

# Examples:
1. insertion
ref: ATG-AACGTA  -> cleaned_ref: ATGAACGTA 
seq: ATGTAACGTA  -> cleaned_seq: ATGAACGTA

2. deletion
ref: ATGAACGTA  -> cleaned_ref: ATGAACGTA
seq: ATG-ACGTA  -> cleaned_seq: ATGNACGTA

**NOTE** We always assume the readingFrame is 1
"""
function clean_frameshifts(aligned_ref::LongDNA{4}, aligned_seq::LongDNA{4}; verbose::Bool=false)
    # get ref sequence
    ref = ungap(aligned_ref)
    # exception handling
    length(aligned_ref) == length(aligned_seq) || throw(ArgumentError("aligned sequences have differing lengths!"))
    length(ref) % 3 == 0 || throw(ArgumentError("The original reference sequence (ungap(aligned_ref)) must have length divisible by 3"))
    !any(
        (base_ref == DNA_Gap) && (base_seq == DNA_Gap) for (base_ref,base_seq) in zip(aligned_ref, aligned_seq)
    ) || throw(ArgumentError("Both sequences have a gap at the same index, which therefore cannot be resolved"))
    # speed_up by checking if any frameshifts are present in the first place
    if !_has_frameshifts(aligned_ref, aligned_seq)
        return aligned_ref, aligned_seq
    else
        return _clean_frameshifts(aligned_ref, aligned_seq, verbose=verbose)
    end
end
# internal _clean_frameshifts method
function _clean_frameshifts(aligned_ref::LongDNA{4}, aligned_seq::LongDNA{4}; verbose::Bool=false)
    # codon bookkeeping variables
    insertAddon = 0
    ref = ungap(aligned_ref)
    codon_length = length(ref)÷3
    # only used in verbose mode
    total_num_of_clean_up_operations = 0
    total_num_deletion_gaps_cleaned = 0
    total_num_insertion_gaps_cleaned = 0
    # boolean lambda functions
    is_3gap = seq -> (seq == LongDNA{4}("---"))
    no_gaps = seq -> !any(==(DNA_Gap), seq)
    # initialize empty alignment
    cleaned_ref = LongDNA{4}("")
    cleaned_seq = LongDNA{4}("")
    codon_index = 1
    last_alignment_index = 0
    while codon_index <= codon_length
        cur_pos = 3*(codon_index-1)+1+insertAddon
        cur_codon_ref = aligned_ref[cur_pos:cur_pos+2]
        cur_codon_seq = aligned_seq[cur_pos:cur_pos+2]
        if (no_gaps(cur_codon_ref) && no_gaps(cur_codon_seq)) || is_3gap(cur_codon_seq)
            append!(cleaned_ref, cur_codon_ref)
            append!(cleaned_seq, cur_codon_seq)
            codon_index += 1
            last_alignment_index = cur_pos+2
            continue
        elseif is_3gap(cur_codon_ref)
            insertAddon += 3
            append!(cleaned_ref, cur_codon_ref)
            append!(cleaned_seq, cur_codon_seq)
            last_alignment_index = cur_pos+2
            continue
        end
        # for visibility when - verbose
        num_insertion_gaps_cleaned = 0
        num_deletion_gaps_cleaned = 0
        # look through codon and try to fix frameshift
        for j in 1:3
            # insertions reletive to reference are resolved by removal
            if cur_codon_ref[j] == DNA_Gap
                num_insertion_gaps_cleaned += 1
            # deletions relative to reference are resolved by inserting ambigious nucleotide N
            elseif cur_codon_seq[j] == DNA_Gap
                num_deletion_gaps_cleaned += 1
                push!(cleaned_ref, cur_codon_ref[j])
                append!(cleaned_seq, LongDNA{4}("N")) # ambigious nucleotide
            else
                push!(cleaned_ref, cur_codon_ref[j])
                push!(cleaned_seq, cur_codon_seq[j])
            end
        end

        # add extra bases to the codon to account for the insertions that were removed
        extra_bases_needed = num_insertion_gaps_cleaned
        new_needed_insertAddon = 0
        num_codon_bases_fixed = 0
        while num_codon_bases_fixed < extra_bases_needed
            new_needed_insertAddon += 1
            if aligned_ref[cur_pos+2+new_needed_insertAddon] == DNA_Gap
                num_insertion_gaps_cleaned += 1
                continue
            elseif aligned_seq[cur_pos+2+new_needed_insertAddon] == DNA_Gap
                num_deletion_gaps_cleaned += 1
                num_codon_bases_fixed += 1
                push!(cleaned_ref, aligned_ref[cur_pos+2+new_needed_insertAddon])
                append!(cleaned_seq, LongDNA{4}("N")) # ambigious nucleotide
            else
                num_codon_bases_fixed += 1
                push!(cleaned_ref, aligned_ref[cur_pos+2+new_needed_insertAddon])
                push!(cleaned_seq, aligned_seq[cur_pos+2+new_needed_insertAddon])
            end
        end
        last_alignment_index = cur_pos+2+new_needed_insertAddon
        # if we get 3 ambigious nucleotides inserted in a codon we replace it with the original codon deletion
        if cleaned_seq[length(cleaned_ref)-2:end] == LongDNA{4}("NNN")
            num_deletion_gaps_cleaned -= 3
            cleaned_seq = cleaned_seq[1:length(cleaned_ref)-3] * LongDNA{4}("---")
        end

        # verbose visibility of edits made to the alignment
        if verbose && (num_insertion_gaps_cleaned != 0 || num_deletion_gaps_cleaned != 0)
            # update operation statistics
            total_num_deletion_gaps_cleaned += num_deletion_gaps_cleaned
            total_num_insertion_gaps_cleaned += num_insertion_gaps_cleaned
            # total number of operations across clean_up
            total_num_of_clean_up_operations = total_num_deletion_gaps_cleaned + total_num_insertion_gaps_cleaned
            # visibility for changes at current codon
            if num_insertion_gaps_cleaned != 0
                println("removed $num_insertion_gaps_cleaned insertion(s) at codon position $codon_index relative to reference")
            end
            if num_deletion_gaps_cleaned != 0
                println("inserted $num_deletion_gaps_cleaned ambigious nucleotide(s) N at codon position $codon_index relative to reference")
            end
            if codon_index == 1
                println("left codon is edited!")
                println("ref: ", ref[1:3*(codon_index+1)])
                println("original alignment_ref: ", aligned_ref[cur_pos:last_alignment_index])
                println("original alignment_seq: ", aligned_seq[cur_pos:last_alignment_index])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-2:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-2:end])
                println("_____________________________________________________")
            elseif codon_index == codon_length
                println("right codon is edited!")
                println("ref: ", ref[3*(codon_index-2)+1:3*(codon_index)])
                println("original alignment_ref: ", aligned_ref[cur_pos-3:last_alignment_index])
                println("original alignment_seq: ", aligned_seq[cur_pos-3:last_alignment_index])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-5:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-5:end])
                println("_____________________________________________________")
            else
                println("middle codon is edited!")
                println("ref: ", ref[3*(codon_index-2)+1:3*(codon_index+1)])
                println("original alignment_ref: ", aligned_ref[cur_pos-3:last_alignment_index])
                println("original alignment_seq: ", aligned_seq[cur_pos-3:last_alignment_index])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-5:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-5:end])
                println("_____________________________________________________")
            end
        end
        insertAddon += new_needed_insertAddon
        codon_index += 1
    end
    # check if insertion occurs after last codon in reference sequence
    num_insertion_gap_after_last_codon = length(aligned_ref) - last_alignment_index
    if num_insertion_gap_after_last_codon != 0
        num_extra_insertion_gaps_cleaned = num_insertion_gap_after_last_codon % 3
        append!(cleaned_ref, aligned_ref[last_alignment_index+1:length(aligned_ref)-num_extra_insertion_gaps_cleaned])
        append!(cleaned_seq, aligned_seq[last_alignment_index+1:length(aligned_ref)-num_extra_insertion_gaps_cleaned])
        if verbose && (num_extra_insertion_gaps_cleaned != 0)
            # add the final operations to the tally
            total_num_insertion_gaps_cleaned += num_extra_insertion_gaps_cleaned
            total_num_of_clean_up_operations += num_extra_insertion_gaps_cleaned
            # visibility printout
            println("the last $num_extra_insertion_gaps_cleaned insertions at the end of the alignment have been removed")
        end
    end
    # explain results of clean_up
    if verbose
        if total_num_of_clean_up_operations != 0
            println("Alignment had frameshift mutations:\nIn total $total_num_insertion_gaps_cleaned insertion(s) were removed and $total_num_deletion_gaps_cleaned ambigious nucleotide(s) were added to non-reference sequence")
        else
            println("Alignment has no frameshift mutations - no clean up needed.")
        end
    end
    return cleaned_ref, cleaned_seq
end
"""
clean a multiple sequence alignment provided one of the sequence is a reference sequence
"""
function clean_frameshifts(aligned_ref::LongDNA{4}, aligned_seqs::Vector{LongDNA{4}})
    # exception handling
    !any(
        length(aligned_ref) != length(seq) for seq in aligned_seqs
    ) || throw(ArgumentError("Invalid input:\n Not all of the aligned sequences in given MSA have the same length!"))
    # clean_up heuristic - clean pairwise and scaffold
    N = length(aligned_seqs)
    cleaned_refs = Vector{LongDNA{4}}(undef, N)
    cleaned_seqs = Vector{LongDNA{4}}(undef, N)
    for i in 1:N
        # project msa to valid pairwise alignment by removing redundant gaps - i.e. gaps matched to gaps
        stripped_aligned_ref, stripped_aligned_seq = strip_gap_only_cols(aligned_ref, aligned_seqs[i])
        # then clean the new pairwise alignment
        cleaned_ref, cleaned_seq = clean_frameshifts(stripped_aligned_ref, stripped_aligned_seq)
        # save result
        cleaned_refs[i] = cleaned_ref
        cleaned_seqs[i] = cleaned_seq
    end
    # scaffold the cleaned pairwise alignments into new frameshift-free multiple sequence alignment
    cleaned_msa = scaffold_msa_from_pairwise(cleaned_refs, cleaned_seqs) # TODO poorly optimized
    return cleaned_msa
end

function strip_gap_only_cols(seq1::LongDNA{4}, seq2::LongDNA{4})
    len = length(seq1)
    # count how many columns should be removed
    cnt = 0
    @inbounds for i in 1:len
        if !(seq1[i] == DNA_Gap && seq2[i] == DNA_Gap)
            cnt += 1
        end
    end
    # return original output if no cols contains only gaps
    if cnt == 0
        return seq1, seq2
    else
        # preallocate container
        out1 = LongDNA{4}("-"^cnt)
        out2 = LongDNA{4}("-"^cnt)
        # fill output vectors
        pos = 1
        @inbounds for i in 1:len
            if !(seq1[i] == DNA_Gap && seq2[i] == DNA_Gap)
                out1[pos] = seq1[i]
                out2[pos] = seq2[i]
                pos += 1
            end
        end
        # return output
        return out1, out2
    end
end

@inline function _has_frameshifts(aligned_ref::LongDNA{4}, aligned_seq::LongDNA{4})
    # assumes valid nucleotide alignment with CDS ref
    gap_run = 0
    @inbounds for i in eachindex(aligned_ref)
        if aligned_ref[i] == DNA_Gap || aligned_seq[i] == DNA_Gap
            gap_run += 1
            if gap_run % 3 != 0
                return true
            end
        else
            gap_run = 0
        end
    end
    return false
end