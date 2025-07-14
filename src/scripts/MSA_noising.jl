import Pkg
Pkg.add("Distributions")

using BioSequences
using Distributions

# TODO make this more realistic
function mutateSequence(seq::LongDNA{4}, seqId::String="missing")
    n = length(seq)
    dna = LongDNA{4}("ACGT")
    num_mutations = rand(Poisson(10))
    num_indels = rand(Poisson(10))
    # substritution mutations
    for i in 1:num_mutations
        mutation_pos = rand(2:n-1)
        seq = seq[1:mutation_pos-1] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[mutation_pos+1:end]
    end
    # insertion/deletion mutations
    for i in 1:num_indels
        indel_pos = rand(1:1+(n-4)÷3) # idk why like this
        a = rand(1:10)
        println("cur length: ", n, " ", length(seq))
        println("worst indel pos: ", indel_pos+3)
        if a >= 8
            # insertion
            println("aqui")
            seq = seq[1:3*indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 3) * seq[3*indel_pos + 1 : end]
            println(seqId," indel_pos: ", indel_pos)
            n += 3
        elseif a >= 5
            println("aqui 2")
            seq = seq[1:3*indel_pos] * seq[3*indel_pos + 4 : end]
            println(seqId," indel_pos: ", indel_pos)
            n -= 3
        elseif a > 3
            seq = seq[1:3*indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[3*indel_pos + 1 : end]
            n += 1
        else
            # deletion
            deleteat!(seq, 3*indel_pos+1)
            n -= 1
        end
    end

    return seq
end

function mutateSequence(seq::LongDNA{4}, indel_pos::Int64)
    seq = seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[indel_pos + 1 : end]
end

#    for pos in keys(indel_dict)
#        a, b = indel_dict(pos)
#        # insertion
#        if a == 1
#             
#        end
#
#        if a == -1

#        end
#    end
#end

#function fasta_mutate_sequences(fasta_file::String)
#    seqs, seq_names = fasta_to_LongDNA(fasta_file)
#    # noise all but reference
#    noised_seqs = mutateSequence.(seqs[2:end])
#    write_LongDNA_alignment_to_fasta(vcat(seqs[1],noised_seqs), seq_names, "./fasta_output/noised.fasta")
#end

#function levenshtein_dist
#include("MSA_noising.jl")
#include("seed_chain_align.jl") 
#include("needleman_wunsch.jl")

#collect sequnces and ungap sequence 
#names, refAndseqs = NextGenSeqUtils.read_fasta_with_names("src/P018_subset.fasta")
#numOfSequence = length(refAndseqs)-1
#numOfSequence = 5

#nameArray = [String(names[i]) for i in 1:length(names)]
#println(length(nameArray))
#ref = LongDNA{4}(refAndseqs[1])

#AA_seqs = Array{LongAA}(undef,numOfSequence)
#seqs = Array{LongDNA{4}}(undef,numOfSequence)
#for iseq in 1:numOfSequence
#   seqs[iseq] = LongDNA{4}(refAndseqs[iseq+1])
#   AA_seqs[iseq] = translate(ungap(LongDNA{4}(refAndseqs[iseq+1])))
#end

# add some noise into sequences (not reference)
#for i in 1:length(seqs)-1
#    seqs[i] = mutateSequence(seqs[i])
#end

#=
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

#using Plots
#plotlyjs()  
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
=#