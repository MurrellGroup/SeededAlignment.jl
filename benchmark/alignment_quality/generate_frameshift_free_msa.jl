# add a set seed for the generation of dataset for alignment quality assesment
using Random
Random.seed!(42)

# ATTENTION: Please run from top directory of package!
include("./../noising.jl")
# sequence_length
seqlength = 2001
num_seqs = 50
# Initialize dataset
dataset = Vector{LongDNA{4}}(undef, num_seqs+1)
# ref sequence
ref = generate_random_ref(seqlength÷3)
dataset[1] = ref
# generate shared mutations 
# TODO stack mutations for more complex alignment dataset and introduce some independence
shared_mutation = mutateSequence(ref,
    codon_indel_avg=3.0,
    frameshift_indel_avg=0.0,
    sub_mutation_avg = 0.0
)
# add frameshift errors when delegating to contigs
for j in 1:num_seqs
    dataset[j+1] = copy(shared_mutation)
end
write_fasta("./benchmark/alignment_quality/frameshift_free_msa_raw.fasta", dataset)
msa = msa_codon_align(ref, dataset[2:end])
write_fasta("./benchmark/alignment_quality/msa_codon_align.fasta", msa)
# to compare with nw_align
#msa = msa_codon_align(ref, dataset[2:end], use_seeded=false)
#write_fasta("./benchmark/alignment_quality/msa_codon_align.fasta", msa)