

"""
    msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)

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
function msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)
    cleaned_codon_alignment = Vector{LongDNA{4}}(undef, length(seqs)+1)
    # perform pairwise SeededAlignment for each sequence
    aligned_seqs, aligned_refs = align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)
    # clean indels which violate the reference readingFrame
    cleaned_codon_alignment[1] = ref
    cleaned_codon_alignment[2:end] = clean_alignment_readingframe.(aligned_refs,aligned_seqs)

    return cleaned_codon_alignment
end


"""
    align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)

    Returns a seeded alignment with respect to the reference for each sequence. 

"""

function align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)
    aligned_seqs = Vector{LongDNA{4}}(undef,length(seqs))
    aligned_refs = Vector{LongDNA{4}}(undef,length(seqs))
    # perform seeded alignment for each sequence w.r.t. reference sequence
    for seqId in 1:length(seqs)
        aligned_ref, aligned_seq = seed_chain_align(ref,ungap(seqs[seqId]),moveset,score_params)
        # save entire alignment to clean up later
        aligned_seqs[seqId] = aligned_seq
        aligned_refs[seqId] = aligned_ref
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
            x = x-insertAddon
            if x == 1
                aligned_seq = aligned_seq[x+1:end]
                insertAddon += 1
            elseif x == length(aligned_seq)-insertAddon
                aligned_seq = aligned_seq[1:x-1]
                insertAddon += 1
            else
                aligned_seq = aligned_seq[1:x-1] * aligned_seq[x+1:end]
                insertAddon += 1
            end
        end
    end

    return aligned_seq
end

