"""
    msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::Moveset, score_params::ScoringScheme; use_seeded::Bool=true)

Produces a fast codon aligment with consistent readingFrame. 

# Arguments
- `ref::LongDNA{4}`: The reference DNA sequence, assumed to have an readingFrame which starts at its first nucleotide.
- `seqs::Vector{LongDNA{4}}`: A vector of DNA sequences to align against the reference.
- `moveset::Moveset`: Defines the allowable alignment operations (e.g. match nucleotide, single indel, readingFrame respecting triple indel). Assumes moveset is choosen so some moves respect the readingFrame.
- `score_params::ScoringScheme`: Scoring parameters for codon-aware alignment, including match/mismatch and gap penalties.
    
# Returns
- `Vector{LongDNA{4}}`: A vector of aligned sequences (including the reference), with codon-aware gaps that preserve reading frame consistency.

# Description
Produces a fast multiple codon alignment based on a reference with known readingFrame, a moveset and a ScoringScheme. 
This is done by computing a pairwise alignment with respect to the reference for each sequence and then cleaning up single indel noise.

# Example
```julia
ref = LongDNA{4}("ATGACGTGA")  # Reference with known reading frame
seqs = [LongDNA{4}("ATGTCGTGA"), LongDNA{4}("ATGACGAGA")]
moveset = std_codon_moveset()
score_params = std_scoring()

alignment = msa_codon_align(ref, seqs, moveset, score_params)
```
NOTE: The sequences (ungapped) might not be preserved by the alignment by the cleanup step. Some nucleotides might be inserted or removed.  
"""
function msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}; moveset::Moveset=std_codon_moveset(), scoring::ScoringScheme=std_scoring(), use_seeded=true::Bool, match_codons=true::Bool)
    cleaned_codon_alignment = Vector{LongDNA{4}}(undef, length(seqs)+1)
    # perform pairwise seeded alignment for each sequence and clean indels which violate the reference readingFrame
    cleaned_refs, cleaned_seqs = align_all_to_reference(ref, seqs, moveset, scoring, use_seeded = use_seeded, match_codons = match_codons, do_clean_frameshifts=true)
    # resolve codon_insertion ambigiouity via left-stacking relative to reference
    msa = resolve_codon_insertions(ref, cleaned_refs, cleaned_seqs)
    return msa
end

"""
    align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::Moveset, score_params::ScoringScheme; use_seeded::Bool=true)

    Returns a seeded alignment with respect to the reference for each sequence. 

"""

function align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::Moveset, score_params::ScoringScheme; use_seeded=true::Bool, match_codons=true::Bool,do_clean_frameshifts=true::Bool)
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
        aligned_ref, aligned_seq = align(ref = ref, query = ungap(seqs[seqId]),moveset=moveset,scoring=score_params, match_codons=match_codons, do_clean_frameshifts=do_clean_frameshifts) 
        # save entire alignment to clean up later
        aligned_seqs[seqId] = copy(aligned_seq)
        aligned_refs[seqId] = copy(aligned_ref)
    end
    return aligned_refs, aligned_seqs
end

# FIXME assume no single indels and codon_respecting insertions 0 mod 3
function find_triplet_insertions_codonindex(aligned_ref::LongDNA{4})
    # TODO check for non-codon-respecting indels
    #@show length(aligned_ref)
    #@show length(ungap(aligned_ref))
    #@show aligned_ref
    #@assert length(ungap(aligned_ref)) % 3 == 0
    

    # NOTE codon_pos denotes which codon is in front of the insertion
    # find where the insertions happened relativ to codons in the reference
    codon_length = length(ungap(aligned_ref))÷3 # floor division
    insertion_codonindex_list = zeros(Int64, codon_length+1)
    insertionFlag = false
    cur_insertion_length = 0
    pos = 1
    codon_pos = 1 
    # loop through aligned_ref to find codon insertions
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