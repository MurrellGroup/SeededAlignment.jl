"""
    msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme; use_seeded::Bool=true)

Produces a fast codon aligment with consistent readingFrame. 

# Arguments
- `ref::LongDNA{4}`: The reference DNA sequence, assumed to have an readingFrame which starts at its first nucleotide.
- `seqs::Vector{LongDNA{4}}`: A vector of DNA sequences to align against the reference.
- `moveset::MoveSet`: Defines the allowable alignment operations (e.g. match nucleotide, single indel, readingFrame respecting triple indel). Assumes moveset is choosen so some moves respect the readingFrame.
- `score_params::ScoreScheme`: Scoring parameters for codon-aware alignment, including match/mismatch and gap penalties.
    
# Returns
- `Vector{LongDNA{4}}`: A vector of aligned sequences (including the reference), with codon-aware gaps that preserve reading frame consistency.

# Description
Produces a fast multiple codon alignment based on a reference with known readingFrame, a moveset and a scoreScheme. 
This is done by computing a pairwise alignment with respect to the reference for each sequence and then cleaning up single indel noise.

# Example
```julia
ref = LongDNA{4}("ATGACGTGA")  # Reference with known reading frame
seqs = [LongDNA{4}("ATGTCGTGA"), LongDNA{4}("ATGACGAGA")]
moveset = std_codon_moveset()
score_params = std_codon_scoring()

alignment = msa_codon_align(ref, seqs, moveset, score_params)
```
NOTE: The sequences (ungapped) might not be preserved by the alignment by the cleanup step. Some nucleotides might be inserted or removed.  
"""
function msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme; use_seeded=true)
    cleaned_codon_alignment = Vector{LongDNA{4}}(undef, length(seqs)+1)
    # perform pairwise seeded alignment for each sequence
    aligned_refs, aligned_seqs = align_all_to_reference(ref, seqs, moveset, score_params, use_seeded = use_seeded)
    # clean indels which violate the reference readingFrame
    cleaned_refs, cleaned_seqs = clean_alignment_readingframe(aligned_refs, aligned_seqs)
    # resolve codon_insertion ambigiouity via left-stacking relative to reference
    msa = resolve_codon_insertions(ref, cleaned_refs, cleaned_seqs)
    return msa
end

"""
    align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme; use_seeded::Bool=true)

    Returns a seeded alignment with respect to the reference for each sequence. 

"""

function align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme; use_seeded::Bool = true)
    aligned_seqs = Vector{LongDNA{4}}(undef,length(seqs))
    aligned_refs = Vector{LongDNA{4}}(undef,length(seqs))
    # choose alignment method
    if use_seeded
        align = seed_chain_align
    else
        align = nw_align
    end
    # perform alignment for each sequence w.r.t. reference sequence
    for seqId in 1:length(seqs)
        aligned_ref, aligned_seq = align(ref,ungap(seqs[seqId]),moveset,score_params) 
        # save entire alignment to clean up later
        aligned_seqs[seqId] = copy(aligned_seq)
        aligned_refs[seqId] = copy(aligned_ref)
    end
    return aligned_refs, aligned_seqs
end

"""
    clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4})

    Takes a pairwise alignment of a reference (with known reading frame) and a sequence, and removes single indels which
    don't respect the reference's reading frame. To clean a full codon alignment you need multiple pairwise alignments and can clean up 
    single indels by broadcasting with "." 

    # NOTE We assume the readingFrame is 0 mod 3 with sequences 0 indexed
"""
# TODO make these functions more robust and make a nicer exported version, maybe add clean option to other aligners. 
function clean_alignment_readingframe(aligned_ref::LongDNA{4},aligned_seq::LongDNA{4})
    indelDict = Dict()
    N_codon = LongDNA{4}("AAC")
    # find all single indels, -1 deletion, 1 insertion
    for i in 1:length(aligned_seq)
        # first codon
        if (i <= 3)
            if aligned_seq[i] == DNA_Gap && !(aligned_seq[1] == DNA_Gap && aligned_seq[2] == DNA_Gap && aligned_seq[3] == DNA_Gap)
                indelDict[i] = -1
            end
            if aligned_ref[i] == DNA_Gap && !(aligned_ref[1] == DNA_Gap && aligned_ref[2] == DNA_Gap && aligned_ref[3] == DNA_Gap)
                indelDict[i] = 1
            end
        #last codon
        elseif (i >= length(aligned_seq)-2)
            if aligned_seq[i] == DNA_Gap && !(aligned_seq[end-2] == DNA_Gap && aligned_seq[end-1] == DNA_Gap && aligned_seq[end] == DNA_Gap)
                indelDict[i] = -1
            end
            if aligned_ref[i] == DNA_Gap && !(aligned_ref[end-2] == DNA_Gap && aligned_ref[end-1] == DNA_Gap && aligned_ref[end] == DNA_Gap)
                indelDict[i] = 1
            end
        # any other codon
        else
            if (aligned_seq[i] == DNA_Gap) && !(aligned_seq[i+1] == DNA_Gap || aligned_seq[i-1] == DNA_Gap)
                indelDict[i] = -1
            end
            if (aligned_ref[i] == DNA_Gap) && !(aligned_ref[i+1] == DNA_Gap || aligned_ref[i-1] == DNA_Gap)
                indelDict[i] = 1
            end
        end
    end
    indelIndicies = collect(keys(indelDict))
    indelIndicies_sorted = sort(indelIndicies)

    insertAddon = 0
    for x in indelIndicies_sorted
        # deletion, 
        if indelDict[x] == -1
            x = x-insertAddon
            if x <= 3
                r = (x-1)%3
                startInsertPos = x-r
                aligned_seq = N_codon * aligned_seq[startInsertPos+3:end]
            elseif x >= length(aligned_seq)-2-insertAddon
                r = (x-1)%3
                startInsertPos = x-r
                aligned_seq = aligned_seq[1:startInsertPos-1] * N_codon
            else
                r = (x-1)%3
                startInsertPos = x-r
                aligned_seq = aligned_seq[1:startInsertPos-1] * N_codon * aligned_seq[startInsertPos+3:end]
            end
        # insertion
        else
            #println("insertion removal")
            x = x-insertAddon
            if x == 1
                aligned_seq = aligned_seq[x+1:end]
                aligned_ref = aligned_ref[x+1:end]
                insertAddon += 1
            elseif x == length(aligned_seq)-insertAddon
                aligned_seq = aligned_seq[1:x-1]
                aligned_ref = aligned_ref[1:x-1]
                insertAddon += 1
            else
                #println("insertion removal middle removed: ", aligned_seq[x], " at ", x)
                aligned_seq = aligned_seq[1:x-1] * aligned_seq[x+1:end]
                aligned_ref = aligned_ref[1:x-1] * aligned_ref[x+1:end]
                insertAddon += 1
            end
        end
    end
    return aligned_ref, aligned_seq
end

function clean_alignment_readingframe(aligned_ref_vec::Vector{LongDNA{4}},aligned_seq_vec::Vector{LongDNA{4}})
    @assert length(aligned_ref_vec) == length(aligned_seq_vec)
    n = length(aligned_ref_vec)

    cleaned_refs = Vector{LongDNA{4}}(undef, length(aligned_ref_vec))
    cleaned_seqs = Vector{LongDNA{4}}(undef, length(aligned_seq_vec))
    for i in 1:n
        cleaned_refs[i], cleaned_seqs[i] = clean_alignment_readingframe(aligned_ref_vec[i],aligned_seq_vec[i])
    end
    return cleaned_refs, cleaned_seqs
end

# FIXME assume no single indels and codon_respecting insertions 0 mod 3
# FIXME make sure ref has length divisible by 3
function find_triplet_insertions_codonindex(aligned_ref::LongDNA{4})
    codon_length = length(ungap(aligned_ref))÷3 # floor division
    @assert length(ungap(aligned_ref)) % 3 == 0
    # TODO check for non-codon-respecting indels
    insertion_codonindex_list = zeros(Int64, codon_length+1)
    insertionFlag = false
    cur_insertion_length = 0
    pos = 1
    codon_pos = 1 
    # NOTE codon_pos denotes which codon is in front of the insertion
    # find where the insertions happened relativ to codons in the reference
    while pos <= length(aligned_ref)
        if aligned_ref[pos] == DNA_Gap
            cur_insertion_length += 3
            if !insertionFlag
                insertionFlag = true
            end
        else 
            if insertionFlag
                insertionFlag = false
                insertion_codonindex_list[codon_pos] = cur_insertion_length
                cur_insertion_length = 0
                codon_pos += 1
            else
                codon_pos += 1
            end
        end
        pos += 3
    end
    if cur_insertion_length > 0
        insertion_codonindex_list[end] = cur_insertion_length
    end
    return insertion_codonindex_list
end

# TODO implement the speedup version of this
function resolve_codon_insertions(ref::LongDNA{4}, cleaned_refs::Vector{LongDNA{4}}, cleaned_seqs::Vector{LongDNA{4}})
    # bookkeeping variables
    n_seqs = length(cleaned_seqs)
    ref_insert_addon = 0
    codon_length = (length(ungap(ref))÷3)
    cleaned_codon_alignment = Vector{LongDNA{4}}(undef, n_seqs+1)
    insertion_bookeeping_matrix = zeros(Int64,n_seqs,codon_length+1)
    
    for i in 1:n_seqs
        # backbone of our msa
        cleaned_codon_alignment[i+1] = cleaned_seqs[i]
        # keep track where the codon insertions are happening
        insertion_bookeeping_matrix[i,:] = find_triplet_insertions_codonindex(cleaned_refs[i])
    end
    cleaned_codon_alignment[1] = ref
    # visually fix the codon insertions by that we left-stack codon insertions at each insertion intervall.
    for codon_index in 1:(codon_length+1) # TODO this can be speedup by always concatenating insertions at the end of sequences.
        # check for insertions that end at the codon_index:th codon relativ to reference
        longest_insert = maximum(insertion_bookeeping_matrix[:,codon_index])
        if longest_insert > 0
            if codon_index == 1
                cleaned_codon_alignment[1] = repeat(LongDNA{4}("-"), longest_insert) * cleaned_codon_alignment[1][1:end]
                for seqId in 1:n_seqs
                    insertion_gap_length = insertion_bookeeping_matrix[seqId,codon_index]
                    if !(insertion_gap_length == longest_insert)
                        cleaned_codon_alignment[seqId+1] = (
                        cleaned_codon_alignment[seqId+1][1:insertion_gap_length] * 
                        repeat(LongDNA{4}("-"), longest_insert-insertion_gap_length) * 
                        cleaned_codon_alignment[seqId+1][insertion_gap_length+1:end]
                        )
                    end
                end
                    
            elseif codon_index == codon_length+1
                cleaned_codon_alignment[1] = cleaned_codon_alignment[1][1:end] * repeat(LongDNA{4}("-"), longest_insert)
                for seqId in 1:n_seqs
                    insertion_gap_length = insertion_bookeeping_matrix[seqId,codon_index]
                    if !(insertion_gap_length == longest_insert)
                        cleaned_codon_alignment[seqId+1] = (
                        cleaned_codon_alignment[seqId+1][1:end] * 
                        repeat(LongDNA{4}("-"), longest_insert-insertion_gap_length) 
                        )
                    end
                end

            else
                
                cleaned_codon_alignment[1] = (
                cleaned_codon_alignment[1][1:3*(codon_index-1)+ref_insert_addon] * 
                repeat(LongDNA{4}("-"), longest_insert) *
                cleaned_codon_alignment[1][3*(codon_index-1)+ref_insert_addon+1:end]
                )
                for seqId in 1:n_seqs
                    insertion_gap_length = insertion_bookeeping_matrix[seqId,codon_index]
                    if !(insertion_gap_length == longest_insert)
                        cleaned_codon_alignment[seqId+1] = (
                        cleaned_codon_alignment[seqId+1][1:3*(codon_index-1)+ref_insert_addon+insertion_gap_length] * 
                        repeat(LongDNA{4}("-"), longest_insert-insertion_gap_length) * 
                        cleaned_codon_alignment[seqId+1][1+3*(codon_index-1)+ref_insert_addon+insertion_gap_length:end]
                        )
                    end
                end
            end
            ref_insert_addon += longest_insert
        end
    end

    return cleaned_codon_alignment
end