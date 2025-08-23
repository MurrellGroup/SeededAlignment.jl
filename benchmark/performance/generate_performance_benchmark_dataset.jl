# add a set seed for the generation of performance benchmarking dataset
using Random
Random.seed!(42)

# ATTENTION: Please run from top directory of package!
include("./../noising.jl")
# sequence_length
const seqlength = 2001
const num_seqs = 200
# Initialize dataset
dataset = Vector{LongDNA{4}}(undef, num_seqs+1)
# ref sequence
ref = generate_random_ref(seqlength÷3)
dataset[1] = ref
# generate shared mutations
shared_mutation = mutateSequence(ref,
    codon_indel_avg=7.0,
    frameshift_indel_avg=0.0,
    sub_mutation_avg = 3.0
)
# add frameshift errors when delegating to contigs
for j in 1:num_seqs
    dataset[j+1] = mutateSequence(shared_mutation, codon_indel_avg=1.0, frameshift_indel_avg = 3.0,sub_mutation_avg = 5.0)
end
write_fasta("./benchmark/performance/benchmark_input_sequences.fasta", dataset)
msa = msa_codon_align(ref, dataset[2:end], use_seeded=true)
write_fasta("./benchmark/performance/benchmarked_msa_result.fasta", msa)
# compare with nw_align
#msa = msa_codon_align(ref, dataset[2:end], use_seeded=false)
#write_fasta("./benchmark/performance/benchmarked_msa_result.fasta", msa)