# add a set seed for the generation of performance benchmarking dataset
using Random
Random.seed!(42)

# ATTENTION: Please run from top directory of package!
include("./../noising.jl")

# TODO make it easier to customize similarity and such

# sequence_length
const seqlength = 699
const num_seqs = 30
# Initialize dataset
dataset = Vector{LongDNA{4}}(undef, num_seqs+1)
# ref sequence
ref = generate_random_ref(seqlength÷3)
dataset[1] = ref
for i in 1:num_seqs
    dataset[i+1] = mutateSequence(ref)
end

write_fasta("./benchmark/performance/benchmark_input_sequences.fasta", dataset)