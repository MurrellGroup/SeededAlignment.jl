
"""
    clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4},verbose_flag::Bool=false)

    Takes a pairwise alignment of a reference (with known reading frame) and a sequence, and removes single indels which
    don't respect the reference's reading frame.

    # NOTE We assume the readingFrame is 0 mod 3 with sequences 0 indexed
"""
# TODO add clean option to aligners.
function clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4},verbose_flag::Bool=false)

    @assert length(ungap(aligned_ref)) % 3 == 0 "Ungapped reference sequence length not divisible by 3"
    @assert !any((base_ref == '-') && (base_seq == '-') for (base_ref,base_seq) in zip(aligned_ref, aligned_seq)) "Both sequences have a gap at the same index"

    insertAddon = 0
    codon_length = Int64(length(ungap(aligned_ref))/3)
    ref = ungap(aligned_ref) # only used if verbose_flag
    # boolean lambda functions
    is_3gap = seq -> (seq == LongDNA{4}("---"))
    no_gaps = seq -> (length(ungap(seq)) == length(seq))
    # initialize empty alignment
    cleaned_ref = LongDNA{4}("")
    cleaned_seq = LongDNA{4}("")
    codon_index = 1
    while codon_index <= codon_length
        cur_pos = 3*(codon_index-1)+1+insertAddon
        cur_codon_ref = aligned_ref[cur_pos:cur_pos+2]
        cur_codon_seq = aligned_seq[cur_pos:cur_pos+2]
        if (no_gaps(cur_codon_ref) && no_gaps(cur_codon_seq)) || is_3gap(cur_codon_seq)
            append!(cleaned_ref, cur_codon_ref)
            append!(cleaned_seq, cur_codon_seq)
            codon_index += 1
            continue
        elseif is_3gap(cur_codon_ref)
            insertAddon += 3
            append!(cleaned_ref, cur_codon_ref)
            append!(cleaned_seq, cur_codon_seq)
            continue
        end
        num_type1_cleanup = 0
        num_type2_cleanup = 0
        for j in 1:3
            # insertions reletive to reference are resolved by removal
            if cur_codon_ref[j] == DNA_Gap
                num_type1_cleanup += 1
                continue
            # deletions relative to reference are resolved by inserting ambigious nucleotide N
            elseif cur_codon_seq[j] == DNA_Gap
                num_type2_cleanup += 1
                append!(cleaned_ref, cur_codon_ref[j:j])
                append!(cleaned_seq, LongDNA{4}("N")) # ambigious nucleotide
            else
                append!(cleaned_ref, cur_codon_ref[j:j])
                append!(cleaned_seq, cur_codon_seq[j:j])
            end
        end
        # add extra bases to the codon to account for the insertions that were removed
        extra_bases_needed = num_type1_cleanup
        for j in 1:extra_bases_needed
            append!(cleaned_ref, aligned_ref[cur_pos+2+j:cur_pos+2+j])
            append!(cleaned_seq, aligned_seq[cur_pos+2+j:cur_pos+2+j])
        end
        # verbose visibility of edits made to the alignment
        if verbose_flag && (num_type1_cleanup != 0 || num_type2_cleanup != 0)
            if num_type1_cleanup != 0
                println("removed $num_type1_cleanup insertion(s) at codon position $codon_index relative to reference")
            end
            if num_type2_cleanup != 0
                println("inserted $num_type2_cleanup ambigious nucleotide(s) N at codon position $codon_index relative to reference")
            end
            if codon_index == 1
                println("left codon is edited!")
                println("ref: ", ref[1:3*(codon_index+1)])
                println("original alignment_ref: ", aligned_ref[cur_pos:cur_pos+5])
                println("original alignment_seq: ", aligned_seq[cur_pos:cur_pos+5])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-2:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-2:end])
                println("_____________________________________________________")
            elseif codon_index == codon_length
                println("right codon is edited!")
                println("ref: ", ref[3*(codon_index-2)+1:3*(codon_index)])
                println("original alignment_ref: ", aligned_ref[cur_pos-3:cur_pos+2])
                println("original alignment_seq: ", aligned_seq[cur_pos-3:cur_pos+2])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-5:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-5:end])
                println("_____________________________________________________")
            else
                println("middle codon is edited!")
                println("ref: ", ref[3*(codon_index-2)+1:3*(codon_index+1)])
                println("original alignment_ref: ", aligned_ref[cur_pos-3:cur_pos+5])
                println("original alignment_seq: ", aligned_seq[cur_pos-3:cur_pos+5])
                println("cleaned alignment_ref:  ", cleaned_ref[length(cleaned_ref)-5:end])
                println("cleaned alignment_seq:  ", cleaned_seq[length(cleaned_seq)-5:end])
                println("_____________________________________________________")
            end
        end
        insertAddon += extra_bases_needed
        codon_index += 1
    end
    return cleaned_ref, cleaned_seq
end

function clean_alignment_readingframe(aligned_ref_vec::Vector{LongDNA{4}}, aligned_seq_vec::Vector{LongDNA{4}},verbose_flag::Bool=false)
    @assert length(aligned_ref_vec) == length(aligned_seq_vec)
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