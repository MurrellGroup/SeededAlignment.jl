"""
    msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}; moveset::Moveset=STD_CODON_MOVESET, scoring::ScoringScheme=STD_SCORING, codon_scoring_on=true::Bool)

Computes a visual global MSA (multiple Sequence alignment) of coding sequences based on pairwise alignments to a trusted CDS reference `ref` used as an anchor to 
determine the apprioate reading frame coordinates for the other coding sequences (which may contain frameshift errors). Possible frameshifts errors in the 
pairwise alignments are cleaned up and then scaffolded to create a multiple sequence alignment. This is done by left-stacking codon insertions relative to 
the reference. 

Note that this doesn't qualify as a proper multiple sequence alignment in the traditional sense since the aligned sequences are only scored as being 
aligned against the reference sequnece and not each other. 

Even so, it can still provides a useful visualization or approximation for a protein multiple sequence alignment on a nucleotide level. 

# Extended Help

# Arguments
- `ref::LongDNA{4}`: Anchored trusted CDS which decides the reading frame coordinates in the alignment
- `seqs::Vector{LongDNA{4}}`: Coding Sequences (with possible frameshifts due to e.g. sequencing or annotation errors) which are aligned to `ref` and adopts its reading frame coordinates. 
- `moveset::Moveset=STD_CODON_MOVESET`: Defines allowable alignment moves (e.g., codon insertions/deletions).
- `scoring::ScoringScheme=STD_SCORING`: The scoring scheme used during alignment.
- `codon_scoring_on::Bool=true`: Whether to apply additional scoring on codon-level 

# Returns
- `msa::Vector{LongDNA{4}}`: a frameshift-free multiple sequence alignment

# Example

1. codon alignment with no frameshifts present in inputs 

```julia
ref =  LongDNA{4}("ATGTTTCCCGGGTAA")
seq1 = LongDNA{4}("ATGTTTTTTCCCGGGTAA")
seq2 = LongDNA{4}("ATGATGTTTTTTCCCGGGTAAGGG")
seq3 = LongDNA{4}("ATGCCCGGG")

msa = msa_codon_align(ref, [seq1,seq2,seq3])
#= alignment results:

msa[1] == LongDNA{4}("---ATG---TTTCCCGGGTAA---")
msa[2] == LongDNA{4}("---ATGTTTTTTCCCGGGTAA---")
msa[3] == LongDNA{4}("ATGATGTTTTTTCCCGGGTAAGGG")
msa[4] == LongDNA{4}("---ATG------CCCGGG------")

=#
```

"""
function msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}; 
        moveset::Moveset=STD_CODON_MOVESET, 
        scoring::ScoringScheme=STD_SCORING, 
        codon_scoring_on::Bool = true,
        verbose::Bool = false,
        use_seeded::Bool = true
)
    # perform pairwise seeded alignment for each sequence and clean indels which violate the reference readingFrame
    cleaned_refs, cleaned_seqs = align_all_to_reference(ref, seqs, 
        moveset = moveset, 
        scoring = scoring, 
        codon_scoring_on = codon_scoring_on, 
        do_clean_frameshifts = true,
        verbose = verbose,
        use_seeded = use_seeded
    )
    # resolve codon_insertion ambigiouity via left-stacking relative to reference
    msa = scaffold_msa_from_pairwise(cleaned_refs, cleaned_seqs)
    return msa
end

function align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}; 
        moveset::Moveset, 
        scoring::ScoringScheme, 
        codon_scoring_on=true::Bool, 
        do_clean_frameshifts=true::Bool,
        verbose::Bool = false,
        use_seeded::Bool = true
)
    # Initialize space for alignment
    aligned_seqs = Vector{LongDNA{4}}(undef,length(seqs))
    aligned_refs = Vector{LongDNA{4}}(undef,length(seqs))
    if use_seeded
        align = seed_chain_align
    else
        align = nw_align
    end
    # perform alignment for each sequence w.r.t. reference sequence
    for seqId in 1:length(seqs)
        # cleaned pairwise alignments
        aligned_ref, aligned_seq = align(
            ref = ref, 
            query = seqs[seqId], 
            moveset = moveset, scoring = scoring, 
            codon_scoring_on = codon_scoring_on, 
            do_clean_frameshifts = do_clean_frameshifts,
            verbose = verbose
        )
        # save alignment results
        aligned_refs[seqId], aligned_seqs[seqId] = aligned_ref, aligned_seq
    end
    return aligned_refs, aligned_seqs
end

# NOTE assume no single indels and codon_respecting insertions 0 mod 3
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
function scaffold_msa_from_pairwise(cleaned_refs::Vector{LongDNA{4}}, cleaned_seqs::Vector{LongDNA{4}})
    
    ref = ungap(cleaned_refs[1])
    # NOTE: we assume that all cleaned_refs represent the same sequence
    #!any(ref != ungap(cleaned_refs[i]) for i in 2:length(cleaned_refs)) || throw(ArgumentError(""))
    
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
    for codon_index in 1:(codon_length+1) 
        
        # TODO this can be speedup by avoiding concatenation and using push! or append! instead.
        # TODO ALTERNATIVELY try copyto! still probably slower

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