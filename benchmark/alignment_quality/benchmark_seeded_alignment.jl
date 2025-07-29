# alignment quality benchmarking script for comparing the alignment results of seed_chain_align and nw_align

# ATTENTION: Please run from top directory of package!

#= Benchmark Strategy: 
We want to evaluate 3 things. 

MAIN benchmark: This will always be computed
    1. How well do seed_chain_align and nw_align recover the "true protein sequences" after cleaning frameshifts? (MAIN)

EXTRA benchmarks: These are done with try catch and are sometimes ignored to due unhandled case of comparing alignments of differing lengths.
    2. How well does seed_chain_align recover the "true pairwise alignment" on a nucleotide and protein level? (EXTRA)
    3. How well does nw_align/seed_chain_align + cleanup and scaffold strategy recover the "true protein multiple sequence alignment". (EXTRA)
    Note: This is equivalent to assessing the quality of the MSAs made from msa_codon_align. 

Description: 
    For (1) we have a "trusted" frameshift-free multiple sequence alignment and introuduce single-indel noise into each sequence. 
The sequences are then ungapped and aligned via nw_align and seed_chain_align relative to reference. The alignments are cleaned of frameshifts 
and scaffolded into a frameshift free multiple sequence alignment, whereby each sequence of which we compare with the original sequence (unnoised) 
and compute the levenshtein distance between the AA translations of the two sequences. For manual inspection non-homologous regions are 
displayed in printouts and the computed distances are used to produce histograms for both seed_chain_align and nw_align.
    
    For (2) we treat the alignment from nw_align as the "true alignment" but the might be subject to change. We attempt to compute the % of columns matched between the 
two alignments if possible. 
(FIXME: For our TC_scoring we assume that both alignments have the same lengths which isn't might break the scoring in some cases).
    
    For (3) we use the 3 MSAs from (1) - these being the "true" frameshift-free multiple sequence alignment and the MSAs made from nw_align and
seed_chain_align respectively. We attempt to compute the SP_score and TC_score for each of the alignments. 
(FIXME: For our TC_scoring we assume that both alignments have the same lengths which isn't might break the scoring in some cases).
=#

# include noising methods
include("./../noising.jl")
# alignment scoring
include("./alignment_distances.jl")
# set the seed
using Random
Random.seed!(44)

#= 1. How well do seed_chain_align and nw_align recover the "true protein sequences" after cleaning frameshifts?
________________________________________________________________________________________________________________=#
using Plots

#collect sequnces and ungap sequence 
names, ref_and_seqs = read_fasta(".fasta_input/P018_subset.fasta")
ref = ungap(ref_and_seqs[1])
seqs = [ungap(seq) for seq in ref_and_seqs[2:end]]
num_seqs = length(seqs)
# convert seqs to AA_seqs
AA_seqs = LongAA[translate(seqs[i]) for i in 1:num_seqs]
# noise the sequences bad noise
noised_seqs = LongDNA{4}[frameshift_noise_seq!(seqs[i]) for i in 1:num_seqs]
msa = msa_codon_align(ref, noised_seqs)
write_fasta(".fasta_output/msa_problem.fasta", msa)
noised_AA_seqs = LongAA[translate(ungap(msa[i+1])) for i in 1:num_seqs]
levenshtein_distances = Int64[levenshtein(noised_AA_seqs[i], AA_seqs[i]) for i in 1:num_seqs]
# Create a histogram
hist = histogram(levenshtein_distances, bins=5, xlabel="AA_levenshtein_distance", ylabel="hits", title="seed_chain_align - protein seq denoising")
display(hist)
savefig("seed_histogram.png")

# show differences
for i in 1:num_seqs
    cur_denoised = noised_AA_seqs[i]
    cur_ref = AA_seqs[i]
    @assert length(cur_denoised) == length(cur_ref)
    for j in 1:length(cur_ref)
        if cur_ref[j] != cur_denoised[j] && (cur_denoised[j] != AminoAcid('X'))
            println("seqId: ",i)
            if j-5 > 0 && length(AA_seqs[i]) > j+5
                println("noise_corrected: ",cur_denoised[j-5:j+5])
                println("no noise:        ", cur_ref[j-5:j+5])
            elseif length(AA_seqs[i]) < j+5
                println("noise_corrected: ",cur_denoised[j-5:end])
                println("no noise: ",cur_ref[j-5:end])
            else
                println("noise_corrected: ",cur_denoised[1:j+5])
                println("no noise: ",cur_ref[1:j+5])
            end
        end
    end 
end

#= 2. How well does seed_chain_align recover the "true alignment" on a nucleotide level? (EXTRA)
________________________________________________________________________________________________________________=#

#collect sequnces and ungap sequence 
names, ref_and_seqs = read_fasta(".fasta_input/P018_subset.fasta")
ref = ungap(ref_and_seqs[1])
seq = ungap(ref_and_seqs[end])
target_alignment = nw_align(ref=ref, query=seq)
src_alignment = seed_chain_align(ref=ref, query=seq)
write_fasta(".fasta_output/pairwise_comparison.fasta", (target_alignment..., src_alignment...))
try
    nuc_matched_percentage = SP_score(target_alignment, src_alignment)
    println("(nucleotide-level) alignment column match %: ", nuc_matched_percentage)
    # TODO also do SP_score on protein-level
    # TODO investiate if seeding quality can be improved
catch
    println("alignments of different dimensions couldn't be compared.\nSkipping: 2. How well does seed_chain_align recover the 'true alignment' on a nucleotide level? (EXTRA)")
    println("target length: ", length(target_alignment[1]))
    println("src length: ", length(src_alignment[1]))
end

#=3. How well does nw_align/seed_chain_align + cleanup and scaffold strategy recover the "true protein multiple sequence alignment". (EXTRA)
Note: This is equivalent to assessing the quality of the MSAs made from msa_codon_align.
_________________________________________________________________________________________________________________=#

# TODO