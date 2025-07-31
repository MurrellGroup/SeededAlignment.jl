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
(TODO: For our TC_scoring we assume that both alignments have the same lengths which isn't might break the scoring in rare cases).
    
    For (3) we use the 3 MSAs from (1) - these being the "true" frameshift-free multiple sequence alignment and the MSAs made from nw_align and
seed_chain_align respectively. We attempt to compute the SP_score and TC_score for each of the alignments. 
(TODO: For our TC_scoring we assume that both alignments have the same lengths which isn't might break the scoring in rare cases).
=#
using Plots
using BioSequences
# include noising methods
include("./../noising.jl")
# alignment scoring
include("./alignment_distances.jl")
# set the seed
using Random
Random.seed!(44)
# translate nuc_alignment to protein_alignment function
function translate_aligned_nuc_seq(nuc_seq::LongDNA{4})
    @assert length(nuc_seq) % 3 == 0
    # assumes frameshift-free
    codon_len = length(nuc_seq) ÷ 3
    aa_seq = LongAA("")
    isCodonGap = x -> !any([y != DNA_Gap for y in x])
    for i in 1:codon_len
        cur_codon = nuc_seq[3*(i-1)+1:3*i]
        if isCodonGap(cur_codon)
            push!(aa_seq, AA_Gap)
        else
            push!(aa_seq, SeededAlignment.fast_translate((cur_codon[1],cur_codon[2],cur_codon[3])))
        end
    end
    return aa_seq
end
#= 1. How well do seed_chain_align and nw_align recover the "true protein sequences" after cleaning frameshifts?
________________________________________________________________________________________________________________=#

# TODO use public dataset

#collect sequnces and ungap sequence 
names, ref_and_seqs = read_fasta(".fasta_input/P018_subset.fasta")
ref = ungap(ref_and_seqs[1])
seqs = [ungap(seq) for seq in ref_and_seqs[2:end]]
num_seqs = length(seqs)
# convert seqs to AA_seqs
AA_seqs = LongAA[BioSequences.translate(seqs[i]) for i in 1:num_seqs]
# noise the sequences bad noise
noised_seqs = LongDNA{4}[frameshift_noise_seq!(seqs[i]) for i in 1:num_seqs]
msa = msa_codon_align(ref, noised_seqs)
write_fasta(".fasta_output/msa_problem.fasta", msa)
noised_AA_seqs = LongAA[BioSequences.translate(ungap(msa[i+1])) for i in 1:num_seqs]
levenshtein_distances = Int64[levenshtein(noised_AA_seqs[i], AA_seqs[i]) for i in 1:num_seqs]
# Create a histogram
hist = histogram(levenshtein_distances, bins=5, xlabel="AA_levenshtein_distance", ylabel="hits", title="seed_chain_align - protein seq denoising")
display(hist)
savefig("seed_histogram.png")
# show differences 
# TODO we get bad results sometimes... investigate if this is just because the noise is badly done or something deeper
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
# TODO add seeds only comparison. 
target_alignment = nw_align(ref=ref, query=seq)
src_alignment = seed_chain_align(ref=ref, query=seq)
write_fasta(".fasta_output/pairwise_comparison.fasta", (target_alignment..., src_alignment...))
try
    # translate to protein alignment
    prot_target_alignment = translate_aligned_nuc_seq.(target_alignment)
    prot_src_alignment = translate_aligned_nuc_seq.(src_alignment)
    # get nucleotide and protein SP_score
    nuc_matched_percentage = SP_score(target_alignment, src_alignment)
    prot_matched_percentage = SP_score(prot_target_alignment, prot_src_alignment)
    println("(nucleotide-level) alignment column match %: ", nuc_matched_percentage)
    println("(protein-level) alignment column match %: ", prot_matched_percentage)
    # TODO investiate if seeding quality can be improved
catch
    println("alignments of different dimensions couldn't be compared.\nSkipping: 2. How well does seed_chain_align recover the 'true alignment' on a nucleotide level? (EXTRA)")
    println("target length: ", length(target_alignment[1]))
    println("src length: ", length(src_alignment[1]))
end

#=3. How well does nw_align/seed_chain_align + cleanup and scaffold strategy recover the "true protein multiple sequence alignment". (EXTRA)
Note: This is equivalent to assessing the quality of the MSAs made from msa_codon_align.
_________________________________________________________________________________________________________________=#

# TODO implement SP_score and TC_score functions for msa and msa visability. 
# Note this might only be implement once ParitialOrderAlignment msa is working to compare alignment quality.