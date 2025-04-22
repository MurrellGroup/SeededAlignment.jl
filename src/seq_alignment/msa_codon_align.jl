

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
function msa_codon_align(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme, debug=false)
    cleaned_codon_alignment = Vector{LongDNA{4}}(undef, length(seqs)+1)
    # perform pairwise SeededAlignment for each sequence
    aligned_refs, aligned_seqs = align_all_to_reference(ref::LongDNA{4}, seqs::Vector{LongDNA{4}}, moveset::MoveSet, score_params::ScoreScheme)
    # clean indels which violate the reference readingFrame
    cleaned_codon_alignment[1] = ref
    for i in 1:length(seqs)
        cleaned_codon_alignment[i+1] = clean_alignment_readingframe(aligned_refs[i],aligned_seqs[i])
    end
    # handle readingFrame respecting triplet insertion relative to the readingFrame. 
    if !debug
        insertion_dict = Dict{Int64, Set{Int64}}()
        for seqId in 1:length(seqs)
            println("cur seqId: ", seqId)
            triplet_insertions = find_triplet_insertions(aligned_refs[seqId])
            println("tripplet insertion: ",triplet_insertions)
            if length(triplet_insertions) != 0
                for x in triplet_insertions
                    push!(get!(insertion_dict,x, Set{Int64}()), seqId+1) # +1 ignore reference refactor would be good
                end
            end
        end
        println("insertion_dict ", insertion_dict)
        insert_positions = sort!(collect(keys(insertion_dict)))
        println("keys", insert_positions)
        insert_addons_vec = zeros(Int64, length(seqs)+1)
        #inserted_against = fill(0, row_boundary,column_boundary)
        insertAddon = 0
        for insert_pos in insert_positions
            seqsIds_with_insertions = insertion_dict[insert_pos]
            println(" with_insertion", seqsIds_with_insertions)
            # insert triplet if not in the above set
            for seqId in 1:length(seqs)+1
                if !(seqId in seqsIds_with_insertions)
                    println("insertion gaps added", " seqId: ", seqId)
                    #insertAddon = insert_addons_vec[seqId]
                    println("insertAddon: ", insertAddon)
                    cleaned_codon_alignment[seqId] = cleaned_codon_alignment[seqId][1:insert_pos+insertAddon] * LongDNA{4}("---") * cleaned_codon_alignment[seqId][insert_pos+insertAddon+1:end]
                    
                    #insert_addons_vec[seqId] += 3 # could work
                end
            end
            insertAddon += 3
        end
    end
   
    return cleaned_codon_alignment
end

function find_triplet_insertions(aligned_ref::LongDNA{4})
    indicies = Vector{Int64}()
    num_single_insertion = 0
    for i in 1:length(aligned_ref)-2
        if i == 1
            println("early")
        elseif i >= length(aligned_ref)-2
            println("late")
        else
            if aligned_ref[i] == DNA_Gap && !(aligned_ref[i+1] == DNA_Gap || aligned_ref[i-1] == DNA_Gap)
                num_single_insertion += 1
            elseif aligned_ref[i] == DNA_Gap && aligned_ref[i+1] == DNA_Gap && aligned_ref[i+2] == DNA_Gap
                println("triplet identified!")
                println(aligned_ref[i:i+5], " ", i)
                push!(indicies, i-num_single_insertion-1)
            end
        end
    end
    return indicies 
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
        # TODO empty chains bug fix
        #aligned_ref, aligned_seq = seed_chain_align(ref,ungap(seqs[seqId]),moveset,score_params) 
        aligned_ref, aligned_seq = nw_align(ref,ungap(seqs[seqId]),moveset,score_params)
        # save entire alignment to clean up later
        aligned_seqs[seqId] = copy(aligned_seq)
        aligned_refs[seqId] = copy(aligned_ref)
    end
    write_fasta("fasta_output/test_mutated_refs_only2.fasta", aligned_refs)
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
            println("insertion removal")
            x = x-insertAddon
            if x == 1
                aligned_seq = aligned_seq[x+1:end]
                insertAddon += 1
            elseif x == length(aligned_seq)-insertAddon
                aligned_seq = aligned_seq[1:x-1]
                insertAddon += 1
            else
                println("insertion removal middle removed: ", aligned_seq[x], " at ", x)
                aligned_seq = aligned_seq[1:x-1] * aligned_seq[x+1:end]
                insertAddon += 1
            end
        end
    end
    return aligned_seq
end

