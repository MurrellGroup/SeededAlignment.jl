using BioSequences
using SeededAlignment

_, seqs = read_fasta("./benchmark/performance/benchmark_input_sequences.fasta")
# remove gaps 
seqs = LongDNA{4}[ungap(seq) for seq in seqs]
ref = seqs[1]
msa = msa_codon_align(ref, seqs[2:end])
write_fasta("./benchmark/performance/benchmarked_msa_result.fasta", msa)