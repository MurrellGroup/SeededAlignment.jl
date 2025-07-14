include("../MSA_noising.jl")
names, seqs = read_fasta("fasta_input/P018_subset.fasta")
names = [string(i) for i in 1:length(seqs)]
ref = seqs[1]
"""seq2 = mutateSequence(seqs[35])
write_fasta("fasta_output/test_write.fasta", [seqs[1],seq2], seq_names=names[1:2])
# aligntest

println(moveset)
#alignment = seed_chain_align(seqs[1],mutateSequence(seqs[35]),moveset,score_params)
alignment1 = seed_chain_align(seqs[1],seq2,moveset,score_params)[1]
alignment2 = seed_chain_align(seqs[1],seq2,moveset,score_params)[2]
alignment3 = nw_align(seqs[1],seq2,moveset,score_params)[1]
alignment4 = nw_align(seqs[1],seq2,moveset,score_params)[2]
a = clean_alignment_readingframe(alignment1,alignment2)
alignment = [seqs[1],alignment1, alignment2, alignment3, alignment4, a]
write_fasta("fasta_output/test_write.fasta", alignment)
"""
# method of storing mutations
m1 = Move(step=3, score=2.0, stride=3, phase=0, extensionAble=true)
m3 = Move(step=3, score=2.0, stride=3, phase=0, extensionAble=true)
m2 = Move(step=1, score=4.0, stride=1, phase=0, extensionAble=false) # making extensionable can cause problems for some reason...
score_params = ScoreScheme(match_score=0.0, mismatch_score=0.5, extension_score=0.7, edge_ext_begin=true, edge_ext_end=true, kmerlength=12, codon_match_bonus=-2.0)
moveset = MoveSet(match_moves=[Move(3,0.0),Move(1,0.0)], hor_moves=[m3,m2], vert_moves=[m1,m2])

#mutated = mutateSequence.(seqs[2:end], [string(i) for i in 2:length(seqs)])
#write_fasta("fasta_input/mutated_test_seqs.fasta", mutated)

_, mutated = read_fasta("fasta_input/mutated_test_seqs.fasta")

# FIXME next time continue with same mutations via using test_pairwise_fasta as input unmuted. Don't change anything
# we want to fix the issue where we have 2 tripplet insertions at the "same positions"

# FIXME seed_chain_align can fail if no single indel move. 
# TODO get better behaviour when seeding...

#TODO when done an you exten the reference you need to move your deletions relative to reference by 3 since they are offset

#aligned_refs, aligned_seqs = SeededAlignment.align_all_to_reference(ref, 
                                                                    #mutated, 
                                                                    #moveset, 
                                                                    #score_params, 
                                                                    #use_seeded = false, 
                                                                    #codon_matching_enabled=true)
#alignment = nw_align(ref,mutated[1],moveset,score_params)

#a = Vector{LongDNA{4}}()
#b = []
#for i in 1:length(mutated)
#    push!(a, aligned_refs[i])
#    push!(a, aligned_seqs[i])
#end

#for i in 1:length(mutated)
#    push!(b, "ref "*string(i))
#    push!(b, "seq "*string(i))
#end

#display(alignment[3][2015:2030,2015:2030])
#println(alignment[3][end,end])
#aa = clean_alignment_readingframe(alignment[1],alignment[2])

#a = seed_chain_align(ref, mutated[1], moveset, score_params, clean_up_flag=true, codon_matching_enabled=true)
#a = nw_align(ref, mutated[1], moveset, score_params, clean_up_flag=true, codon_matching_enabled=false)

#using BenchmarkTools
ref = ref
mutated = mutated
#a1 = @benchmark nw_align(ref, mutated, moveset, score_params, clean_up_flag=false, codon_matching_enabled=false)
#a2 = @benchmark nw_align(ref, mutated, moveset, score_params, clean_up_flag=true, codon_matching_enabled=false)
#a3 = @benchmark nw_align(ref, mutated, moveset, score_params, clean_up_flag=false, codon_matching_enabled=true)
#display(a1)
#display(a2)
#display(a3)
a = msa_codon_align(ref, mutated, moveset, score_params, codon_matching_enabled=true, use_seeded=true)
write_fasta("fasta_output/test_pairwise.fasta", a) 

#seq_names = ["ref "*string(i), "seq "*string(i), for i in 1:1])
#_, mutated = read_fasta("fasta_input/mutated_test_seqs.fasta")
#alignment = msa_codon_align(ref, mutated, moveset, score_params, use_seeded=false) # change to example
#@show alignment
#write_fasta("fasta_output/test_write3.fasta", alignment, seq_names = [i == 1 ? "ref" : string(i-1) for i in 1:length(alignment)])

# test std_codon_moveset
# msa_codon_align(ref, mutated, std_codon_moveset(), std_codon_scoring())

#14 indel_pos: 1443
#14 indel_pos: 555
#14 indel_pos: 555
#14 indel_pos: 2497
#test_seq = seqs[14]

#test_seq = mutateSequence(test_seq, 555)
#test_seq = mutateSequence(test_seq, 555)
#test_seq = mutateSequence(test_seq, 555)
#test_seq = mutateSequence(test_seq, 555)

#write_fasta("fasta_output/test_mutated_seqs.fasta", [ref, test_seq])
#alignment2 = nw_align(ref, test_seq, moveset, score_params)
#alignment1 = seed_chain_align(ref, test_seq, moveset, score_params)
#alignment = msa_codon_align(ref, [test_seq], moveset, score_params)

#cleaned = clean_alignment_readingframe(alignment[1],alignment[2])
#write_fasta("fasta_output/test_write3.fasta", vcat([alignment2[1], alignment2[2], alignment1[1], alignment1[2]], alignment))