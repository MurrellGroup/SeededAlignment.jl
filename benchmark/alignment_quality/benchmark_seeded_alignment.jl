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

    For (3) we use the 3 MSAs from (1) - these being the "true" frameshift-free multiple sequence alignment and the MSAs made from nw_align and
seed_chain_align respectively. We attempt to compute the SP_score and TC_score for each of the alignments.
=#

using Plots
using BioSequences
using SeededAlignment
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
    aa_ambigious = AminoAcid('X')
    isCodonGap = x -> !any([y != DNA_Gap for y in x])
    for i in 1:codon_len
        cur_codon = nuc_seq[3*(i-1)+1:3*i]
        if DNA_N in cur_codon
            push!(aa_seq, aa_ambigious)
        elseif isCodonGap(cur_codon)
            push!(aa_seq, AA_Gap)
        else
            push!(aa_seq, SeededAlignment.fast_translate(cur_codon[1],cur_codon[2],cur_codon[3]))
        end
    end
    return aa_seq
end
#= 1. How well do seed_chain_align and nw_align recover the "true protein sequences" after cleaning frameshifts?
________________________________________________________________________________________________________________=#

#collect sequnces and ungap sequence 
names, ref_and_seqs = read_fasta("./benchmark/alignment_quality/frameshift_free_msa_raw.fasta")
ref = ungap(ref_and_seqs[1])
seqs = [ungap(seq) for seq in ref_and_seqs[2:end]]
num_seqs = length(seqs)
# convert seqs to AA_seqs
AA_seqs = LongAA[BioSequences.translate(seqs[i]) for i in 1:num_seqs]
# noise the sequences bad noise
noised_seqs = LongDNA{4}[frameshift_noise_seq(seqs[i], frameshift_indel_avg = 1.0) for i in 1:num_seqs]
msa = msa_codon_align(ref, noised_seqs, codon_scoring_on=true, use_seeded=true)
write_fasta("msa.fasta",msa)
noised_AA_seqs = LongAA[BioSequences.translate(ungap(msa[i+1])) for i in 1:num_seqs]
levenshtein_distances = Int64[levenshtein(noised_AA_seqs[i], AA_seqs[i]) for i in 1:num_seqs]
# Create a histogram
hist = histogram(levenshtein_distances, bins=5, xlabel="AA_levenshtein_distance", ylabel="hits", title="seed_chain_align - protein seq denoising")
display(hist) #may not be displayed correctly if called in repl
savefig("seed_histogram.png")
# show differences 
for i in 1:num_seqs
    cur_denoised = noised_AA_seqs[i]
    cur_ref = AA_seqs[i]
    @assert length(cur_denoised) == length(cur_ref)
    for j in 1:length(cur_ref)
        if cur_ref[j] != cur_denoised[j] #&& (cur_denoised[j] != AminoAcid('X')) for investigating only unexpected errors
            println("seqId: ",i)
            println("codon_index: ", j)
            if j-5 > 0 && length(AA_seqs[i]) > j+5
                print("noise_corrected: ")
                color_diff(cur_denoised[j-5:j+5], cur_ref[j-5:j+5], second=false)
                print("no noise:        ") 
                color_diff(cur_denoised[j-5:j+5], cur_ref[j-5:j+5], first=false)
            elseif length(AA_seqs[i]) < j+5
                print("noise_corrected: ")
                color_diff(cur_denoised[j-5:end], cur_ref[j-5:end], second=false)
                print("no noise:        ") 
                color_diff(cur_denoised[j-5:end], cur_ref[j-5:end], first=false)
            else
                print("noise_corrected: ")
                color_diff(cur_denoised[1:j+5], cur_ref[1:j+5], second=false)
                print("no noise:        ") 
                color_diff(cur_denoised[1:j+5], cur_ref[1:j+5], first=false)
            end
        end
    end 
end

seq = noised_seqs[18]
target_alignment = nw_align(ref=ref, query=seq, do_clean_frameshifts=true)
src_alignment = seed_chain_align(ref=ref, query=seq, do_clean_frameshifts=true)
write_fasta(".fasta_output/pairwise_comparison.fasta", (target_alignment..., src_alignment...))
#= 2. How well does seed_chain_align recover the "true alignment" on a nucleotide/codon level? (EXTRA)
________________________________________________________________________________________________________________=#
#collect sequnces and ungap sequence 
println("the following compares the similarity of the pairwise alignments from nw_align and seed_chain_align")
names, ref_and_seqs = read_fasta("./benchmark/alignment_quality/msa_codon_align.fasta")
num_seqs = length(ref_and_seqs)-1
ref = ungap(ref_and_seqs[1])
list_of_similarity_scores_nuc = Vector{Float64}(undef, num_seqs)
list_of_similarity_scores_prot = Vector{Float64}(undef, num_seqs)
for idx in 1:num_seqs
    println("index: ", idx)
    seq = ungap(ref_and_seqs[idx])
    frameshift_noise_seq!(seq)
    target_alignment = nw_align(ref=ref, query=seq, do_clean_frameshifts=true)
    src_alignment = seed_chain_align(ref=ref, query=seq, do_clean_frameshifts=true)
    #write_fasta(".fasta_output/pairwise_comparison.fasta", (target_alignment..., src_alignment...))
    if length(target_alignment[1]) == length(src_alignment[1]) && length(target_alignment[2]) == length(src_alignment[2])
        # translate to protein alignment
        prot_target_alignment = translate_aligned_nuc_seq.(target_alignment)
        prot_src_alignment = translate_aligned_nuc_seq.(src_alignment)
        # get nucleotide and protein SP_score
        nuc_matched_percentage = SP_score(target_alignment, src_alignment)
        println("(nucleotide-level) alignment column match %: ", 100*nuc_matched_percentage, "\n")
        prot_matched_percentage = SP_score(prot_target_alignment, prot_src_alignment)
        println("(protein-level) alignment column match %: ", 100*prot_matched_percentage, "\n")
        list_of_similarity_scores_nuc[idx] = nuc_matched_percentage
        list_of_similarity_scores_prot[idx] = prot_matched_percentage
    else
        println("skipping - alignment dimensions don't match")
        println("target length: ", length(target_alignment[1]))
        println("src length: ", length(src_alignment[1]))
        #list_of_similarity_scores_nuc[idx] = NaN
        list_of_similarity_scores_prot[idx] = NaN
    end
end
# only get stats for proteins 
num_NaN = count(isnan, list_of_similarity_scores_prot)
num_perfect_matches = count(x -> isapprox(x,1.0), list_of_similarity_scores_prot)
filtered_matches = filter(!isnan, list_of_similarity_scores_prot)
median_percent = median(filtered_matches)
min_match_percent = minimum(filtered_matches)
println("Summary information for protein alignments:")
println("num_NaN: ", num_NaN)
println("num_perfect_matches: ", num_perfect_matches)
println("median_percent: ", median_percent)
println("min_match_percent: ", min_match_percent)
argmin_idx = findfirst(x -> x == min_match_percent, list_of_similarity_scores_prot)
println("argmin index: ", argmin_idx)

#=3. How well does nw_align/seed_chain_align + cleanup and scaffold strategy recover the "true protein multiple sequence alignment". (EXTRA)
Note: This is equivalent to assessing the quality of the MSAs made from msa_codon_align.
_________________________________________________________________________________________________________________=#

# TODO implement SP_score and TC_score functions for msa and msa visability.