# performance benchmark script for core methods of SeededAlignment. 

# ATTENTION: Please run from top directory of package!

#= TIPS
1. Ideally results should be at least 1ms-10ms to reduce the influnce of hardware noise between runs
2. Try runing in a low-noise environment to reduce the variance between runs.
3. Note that the runtime of seed_chain_align is strongly influnced by the similarity of the aligned sequences. 
=#
using BenchmarkTools
using SeededAlignment
using BioSequences
using Random

#= TODO this may or may not be realized in the future...
We use 3 seperate datasets which can be edited in generate_performance_benchmark_dataset
    # core functions
    1. dataset for nw_align and seed_chain_align
    2. dataset for msa_codon_align
    3. dataset for clean_frameshifts - msa and pairwise
=#
# retrieve datasets
# _, seqs_pairwise = read_fasta("./benchmark/performance/pairwise.fasta")
# _, seqs_msa_codon_align = read_fasta("./benchmark/performance/msa_codon_align.fasta")
# _, seqs_clean_frameshifts = read_fasta("./benchmark/performance/clean_frameshifts.fasta")
# ungap the inputs if needed
# ref_pairwise, seq_pairwise = ungap(seqs_pairwise[1]), ungap(seqs_pairwise[2])
# ref_msa, seqs_msa = ungap(seqs_msa_codon_align[1]), ungap(seqs_msa_codon_align[2:end]) # make this better and easier for user. 
# ref_clean_frameshifts, seqs_clean_frameshifts = seqs_clean_frameshifts[1], seqs_clean_frameshifts[2:end]

# NOTE dataset contains no gaps
_, seqs = read_fasta("./benchmark/performance/benchmark_input_sequences.fasta")
# remove gaps 
seqs = LongDNA{4}[ungap(seq) for seq in seqs]
ref = seqs[1]
pairwise_seq = seqs[end]
# set up benchmarkGroup
suite = BenchmarkGroup()
# add benchmarks for pairwise alignment methods
suite["seed_chain_align"] = @benchmark seed_chain_align(
    ref=$(ref), 
    query=$(pairwise_seq),
    moveset=$STD_CODON_MOVESET, 
    scoring=$STD_SCORING, 
    codon_scoring_on=$(true), 
    do_clean_frameshifts=$(false), 
    verbose=$(false)
)
suite["nw_align"] = @benchmark nw_align(
    ref=$(ref),
    query=$(pairwise_seq), 
    moveset=$STD_CODON_MOVESET, 
    scoring=$STD_SCORING, 
    codon_scoring_on=$(true), 
    do_clean_frameshifts=$(false), 
    verbose=$(false)
)
# this dataset is very heavy for scaffold_msa due to the number of triplet insertions
# add benchmarks for msa_codon_align
suite["msa_codon_align"] = @benchmark msa_codon_align($ref, $(seqs[2:end]))
# TODO add better dataset to clean. This one has so many insertinos that msa-cleaning is unrealistic. 
msa = msa_codon_align(ref, (seqs[2:end]))
write_fasta("./benchmark/performance/benchmarked_msa_result.fasta", msa)
# add benchmark for cleaning msa 
suite["clean_frameshifts - msa"] = @benchmark clean_frameshifts($msa[1], $(msa[2:end]))
# show results
println("seed_chain_align")
display(suite["seed_chain_align"])
println("nw_align")
display(suite["nw_align"])
println("msa_codon_align")
display(suite["msa_codon_align"])
println("clean_frameshifts - msa")
display(suite["clean_frameshifts - msa"])
# compare with old results if they exist
if isfile("benchmark_core_results.json")
    master_suite = BenchmarkTools.load("benchmark_core_results.json")[1]
else
    @warn "No previous benchmark results found — skipping comparison"
    master_suite = missing
end

if !ismissing(master_suite)

    println("seed_chain_align comparision")
    b1 = median(suite["seed_chain_align"])
    b2 = median(master_suite["seed_chain_align"])
    display(judge(b1, b2))

    println("nw_align")
    b1 = median(suite["nw_align"])
    b2 = median(master_suite["nw_align"])
    display(judge(b1, b2))

    println("msa_codon_align")
    b1 = median(suite["msa_codon_align"])
    b2 = median(master_suite["msa_codon_align"])
    display(judge(b1, b2))

    println("clean_frameshifts - msa")
    b1 = median(suite["clean_frameshifts - msa"])
    b2 = median(master_suite["clean_frameshifts - msa"])
    display(judge(b1, b2))
end

println("To replace the previous baseline with the new result call: update_baseline()")
function update_baseline()
    update_baseline(suite)
end

function update_baseline(new_baseline::BenchmarkGroup)
    BenchmarkTools.save("benchmark_core_results.json", new_baseline)
end

