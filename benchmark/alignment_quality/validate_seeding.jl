# script for quick validatation of which parts of a seeded alignment corresponds to seeds
# code below is an example - modify to match example in mind
using SeededAlignment
using BioSequences

names, ref_and_seqs = read_fasta(".fasta_input/P018_subset.fasta")
ref = ungap(ref_and_seqs[1])
seq = ungap(ref_and_seqs[end])
alignment0 = nw_align(ref=ref, query=seq, codon_scoring_on=true)
alignment1 = seed_chain_align(ref=ref, query=seq, codon_scoring_on=true)
alignment2 = SeededAlignment._seed_chain_align(ref, seq, seed_debug_mode=true)
@assert ungap(alignment0[1]) == ungap(alignment1[1])
@assert ungap(alignment0[2]) == ungap(alignment1[2])
write_fasta("seed.fasta", (alignment0..., alignment1..., alignment2...))