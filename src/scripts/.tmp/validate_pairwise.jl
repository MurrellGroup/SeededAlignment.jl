using BioSequences

_, seqs = read_fasta("fasta_input/P018_subset.fasta")
ref = seqs[1]
_, mutated = read_fasta("fasta_input/mutated_test_seqs.fasta")
B = ungap(mutated[59])

#alignment = nw_align(ref, B, clean_up_flag=true, codon_matching_enabled=true)
alignment = seed_chain_align(ref=ref, query=B, clean_up_enabled=true, codon_matching_enabled=true, verbose=true)

write_fasta("fasta_output/new.fasta", [alignment[1], alignment[2]])

#a = clean_frameshifts(LongDNA{4}("AAATA------TACC"), 
                  #    LongDNA{4}("AAATACCCCCCTACC"),verbose=true)
#@show a