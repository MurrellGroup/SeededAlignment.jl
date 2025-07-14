
"""
    clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4},verbose_flag::Bool=false)

    Takes a pairwise alignment of a reference (with known reading frame) and a sequence, and removes single indels which
    don't respect the reference's reading frame.

    # NOTE We assume the readingFrame is 0 mod 3 with sequences 0 indexed
"""
# TODO add verbose_clean_up to aligners
function clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4},verbose_flag::Bool=false)
    # exception handling
    length(ungap(aligned_ref)) % 3 == 0 || throw(ArgumentError("The original reference sequence (ungap(aligned_ref)) must have length divisible by 3")) 
    !any((base_ref == '-') && (base_seq == '-') for (base_ref,base_seq) in zip(aligned_ref, aligned_seq)) || throw(ArgumentError("Both sequences have a gap at the same index, which therefore cannot be resolved")) 
    
    # codon bookkeeping variables
    insertAddon = 0
    codon_length = length(ungap(aligned_ref))÷3
    # only used if verbose_flag
    ref = ungap(aligned_ref)
    total_num_of_clean_up_operations = 0
    total_num_deletion_gaps_cleaned = 0
    total_num_insertion_gaps_cleaned = 0
    # boolean lambda functions
    is_3gap = seq -> (seq == LongDNA{4}("---"))
    no_gaps = seq -> (length(ungap(seq)) == length(seq))
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
                append!(cleaned_ref, cur_codon_ref[j:j])
                append!(cleaned_seq, LongDNA{4}("N")) # ambigious nucleotide
            else
                append!(cleaned_ref, cur_codon_ref[j:j])
                append!(cleaned_seq, cur_codon_seq[j:j])
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
                append!(cleaned_ref, aligned_ref[cur_pos+2+new_needed_insertAddon:cur_pos+2+new_needed_insertAddon])
                append!(cleaned_seq, LongDNA{4}("N")) # ambigious nucleotide
            else
                num_codon_bases_fixed += 1
                append!(cleaned_ref, aligned_ref[cur_pos+2+new_needed_insertAddon:cur_pos+2+new_needed_insertAddon])
                append!(cleaned_seq, aligned_seq[cur_pos+2+new_needed_insertAddon:cur_pos+2+new_needed_insertAddon])
            end
        end
        last_alignment_index = cur_pos+2+new_needed_insertAddon
        # if we get 3 ambigious nucleotides inserted in a codon we replace it with the original codon deletion
        if cleaned_seq[length(cleaned_ref)-2:end] == LongDNA{4}("NNN")
            num_deletion_gaps_cleaned -= 3
            cleaned_seq = cleaned_seq[1:length(cleaned_ref)-3] * LongDNA{4}("---")
        end

        # verbose visibility of edits made to the alignment
        if verbose_flag && (num_insertion_gaps_cleaned != 0 || num_deletion_gaps_cleaned != 0)
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
        if verbose_flag && (num_extra_insertion_gaps_cleaned != 0)
            # add the final operations to the tally
            total_num_insertion_gaps_cleaned += num_extra_insertion_gaps_cleaned
            total_num_of_clean_up_operations += num_extra_insertion_gaps_cleaned
            # visibility printout
            println("the last $num_extra_insertion_gaps_cleaned insertions at the end of the alignment have been removed")
        end
    end
    # explain results of clean_up
    if verbose_flag
        if total_num_of_clean_up_operations != 0
            println("Alignment had frameshift mutations - in total $total_num_insertion_gaps_cleaned insertions were removed and $total_num_deletion_gaps_cleaned ambigious nucleotides added to non-reference sequence")
        else
            println("Alignment has no frameshift mutations - no clean up needed.")
        end
    end
    return cleaned_ref, cleaned_seq
end

# TODO depricate
function clean_alignment_readingframe(aligned_ref_vec::Vector{LongDNA{4}}, aligned_seq_vec::Vector{LongDNA{4}},verbose_flag::Bool=false)
    
    length(aligned_ref_vec) == length(aligned_seq_vec) || throw(ArgumentError("an equal number of references and non-references must be provided"))
    n = length(aligned_ref_vec)

    cleaned_refs = Vector{LongDNA{4}}(undef, length(aligned_ref_vec))
    cleaned_seqs = Vector{LongDNA{4}}(undef, length(aligned_seq_vec))
    for i in 1:n
        if verbose_flag
            println("Cleaning the $i:th sequence: subsequent edits may follow:")
        end
        cleaned_refs[i], cleaned_seqs[i] = clean_alignment_readingframe(aligned_ref_vec[i],aligned_seq_vec[i],verbose_flag)
    end
    return cleaned_refs, cleaned_seqs
end