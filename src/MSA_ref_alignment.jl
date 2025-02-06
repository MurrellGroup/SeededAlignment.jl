using NextGenSeqUtils
using BioSequences

include("MSA_noising.jl")
include("seed_chain_align.jl")
include("needleman_wunsch.jl")

# collect sequnces and ungap sequence 
names, refAndseqs = NextGenSeqUtils.read_fasta_with_names("src/P018_subset.fasta")
numOfSequence = length(refAndseqs)-1
#numOfSequence = 5

nameArray = [String(names[i]) for i in 1:length(names)]
#println(length(nameArray))
ref = LongDNA{4}(refAndseqs[1])

AA_seqs = Array{LongAA}(undef,numOfSequence)
seqs = Array{LongDNA{4}}(undef,numOfSequence)
for iseq in 1:numOfSequence
   seqs[iseq] = LongDNA{4}(refAndseqs[iseq+1])
   AA_seqs[iseq] = translate(ungap(LongDNA{4}(refAndseqs[iseq+1])))
end

#NOTE that there is at least one single indel which might screw reading frame if mutate
# add some noise into sequences (not reference)
for i in 1:length(seqs)-1
    seqs[i] = mutateSequence(seqs[i])
end

# align pairwise via seeding
function msa_ref_alignment(ref, seqs)
    alignment = Array{LongDNA{4}}(undef,length(seqs)+1)
    alignment[1] = ref
    for seqId in 1:length(seqs)
        # seed pairwise alignment with codon respecting triplet moves
        match_moves = [Move(1,.0), Move(3,.0)]
        # important example
        hor_moves =  [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0,true)]
        vert_moves = [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0,true)]
        aligned_ref, aligned_seq = seed_chain_align(ref, ungap(seqs[seqId]), 0.0, 0.5, match_moves, vert_moves, hor_moves, 0.98*(2.0/3), 21)

        # fix readingframes
        println(seqId)
        fixed_aligned_seq = fix_alignment_readingframe(aligned_ref,aligned_seq)
        alignment[seqId+1] = fixed_aligned_seq
    end
    return alignment
end

#TODO fix indel issues occuring at the last codon. 
# NOTE We assume the readingFrame is 0 mod 3
function fix_alignment_readingframe(aligned_ref,aligned_seq)
    # find all single gaps and get their index. Once we have that we insert N in order unsure of how to fix...
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

msa_alignment = msa_ref_alignment(ref, seqs)
#convert msa_alignment into strings
msa_alignment_str = Array{String}(undef,length(msa_alignment))
for i in 1:length(msa_alignment)
    msa_alignment_str[i] = string(msa_alignment[i])
end

NextGenSeqUtils.write_fasta("post-noise-alignment.fasta", msa_alignment_str, names=nameArray)

levenshtein_distances = Array{Int64}(undef,numOfSequence)
for i in 1:numOfSequence
    println(i)
    println(string(translate(ungap(msa_alignment[i+1]))))
    levenshtein_distances[i] = NextGenSeqUtils.levenshtein(string(translate(ungap(msa_alignment[i+1]))),string(AA_seqs[i]))
end
display(levenshtein_distances)

using Plots
plotlyjs()  
# Switch to the PlotlyJS backend

# Create a histogram
histogram(levenshtein_distances, bins=5, xlabel="levenshtein_distances", ylabel="hits", title="noise-reduction-histogram")
savefig("histogram2.png")

for i in 1:numOfSequence
    cur_pre_noise = translate(ungap(msa_alignment[i+1]))
    for j in 1:length(AA_seqs[i])
        if cur_pre_noise[j] != AA_seqs[i][j]
            println("seqId: ",i)
            if j-5 > 0 && length(AA_seqs[i]) > j+5
                println("noise_corrected: ",cur_pre_noise[j-5:j+5])
                println("no noise:        ",AA_seqs[i][j-5:j+5])
            elseif length(AA_seqs[i]) < j+5
                println("noise_corrected: ",cur_pre_noise[j-5:end])
                println("no noise: ",AA_seqs[i][j-5:end])
            else
                println("noise_corrected: ",cur_pre_noise[1:j+5])
                println("no noise: ",AA_seqs[i][1:j+5])
            end
        end
    end 
end

