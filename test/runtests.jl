using Test

# test dependencies
using SeededAlignment
using BioSequences
using Random

include("../benchmark/noising.jl")

@testset "SeededAlignment.jl" begin
    @testset "1. noisy alignment nw_affine" begin
        Random.seed!(42)
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1001)
        B = mutateSequence(A, verbose=false)
        aligned_A, aligned_B = nw_align(A, B)
        # 1.1 alignment doesn't alter the underlying sequences
        @testset "1.1 alignment doesn't alter the underlying sequences" begin
            @test ungap(aligned_A) == A && ungap(aligned_B) == B
        end
        # 1.2 aligned_sequences have the same length
        @testset "1.2 aligned_sequences have the same length" begin
	        @test length(aligned_A) == length(aligned_B)
        end
        # 1.3 start and ending extension
        @testset "1.3 start and ending extension" begin
            C = LongDNA{4}("TTTAAAGGG")
            D = LongDNA{4}("AAAGGGCCC")
            # we test if we get different alignments from tweaking the boolean parameters
            score_params_on =  ScoringScheme(extension_score=-0.4, edge_ext_begin=true, edge_ext_end=true)
            score_params_off = ScoringScheme(extension_score=-0.4, edge_ext_begin=false,edge_ext_end=false)
            move_set = Moveset(
                ref_insertions = (Move(ref=false, step_length=3, score=-40.0, extendable=true),),
                ref_deletions  = (Move(ref=false, step_length=3, score=-40.0, extendable=true),)
            )
            # test start and end extension on
            alignment = nw_align(C,D,moveset = move_set,scoring = score_params_on)
            @test alignment[1] == LongDNA{4}("TTTAAAGGG---")
            @test alignment[2] == LongDNA{4}("---AAAGGGCCC")
            
            # test start and end extension off
            alignment = nw_align(C,D,moveset = move_set,scoring =score_params_off)
            @test alignment[1] == LongDNA{4}("TTTAAAGGG")
            @test alignment[2] == LongDNA{4}("AAAGGGCCC")
        end
        # 1.4 type inferrence
        @testset "1.4 type inferrence" begin
            @test typeof(@inferred nw_align(A,B)) == Tuple{LongDNA{4},LongDNA{4}}
        end
    end

    @testset "2. ref alignment nw_affine" begin
        move_set = Moveset(
                ref_insertions = (Move(ref=true, step_length=3, score=-10.0, extendable=true),),
                ref_deletions =  (Move(ref=true, step_length=3, score=-10.0, extendable=true),)
            )
        Random.seed!(42)
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 330)
        B = mutateSequence(A, codon_indel_avg = 2.0, frameshift_indel_avg = 0.0, sub_mutation_avg = 0.0, verbose=false)
        # reference informed alignment
        aligned_A, aligned_B = nw_align(ref = A, query = B, moveset=move_set)
        # 2.1 alignment doesn't alter the underlying sequences
        @testset "2.1 alignment doesn't alter the underlying sequences" begin
            @test ungap(aligned_A) == A && ungap(aligned_B) == B
        end
        # 2.2 aligned_sequences have the same length
        @testset "2.2 aligned_sequences have the same length" begin
	        @test length(aligned_A) == length(aligned_B)
        end
        # 2.3 indels don't break codon readingframe
        @testset "2.3 indels don't break codon readingframe" begin
        # extension at edges turned off
        # reference informed alignment
        aligned_A, aligned_B = nw_align(ref = A, query = B, moveset=move_set)
        # check if any clean up was needed
        @test ungap(aligned_A) == A
        @test ungap(aligned_B) == B
        end 
    end

    @testset "3. noisy alignment seed_chain_align" begin
        # TODO better to use more similar sequences so it finds more kmers
        Random.seed!(42)
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1001)
        B = mutateSequence(A, verbose=false)
        aligned_A, aligned_B = seed_chain_align(A, B)
        # 3.1 alignment doesn't alter the underlying sequences
        @testset "3.1 alignment doesn't alter the underlying sequences" begin
            @test ungap(aligned_A) == A && ungap(aligned_B) == B
        end
        # 3.2 aligned_sequences have the same length
        @testset "3.2 aligned_sequences have the same length" begin
	        @test length(aligned_A) == length(aligned_B)
        end
        # 3.3 type inferrence
        @testset "3.3 type inferrence" begin
            @test typeof(@inferred seed_chain_align(A,B)) == Tuple{LongDNA{4},LongDNA{4}}
        end

    end

    @testset "4. ref alignment seed_chain_align" begin
        move_set = Moveset(
                ref_insertions = (Move(ref=true, step_length=3, score=-10.0, extendable=true),),
                ref_deletions =  (Move(ref=true, step_length=3, score=-10.0, extendable=true),)
            )
        Random.seed!(42) # TODO better to use more similar sequences so it finds more kmers
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 330)
        B = mutateSequence(A, codon_indel_avg = 2.0, frameshift_indel_avg = 0.0, sub_mutation_avg = 0.0, verbose=false)
        # reference informed alignment
        aligned_A, aligned_B = seed_chain_align(ref = A, query = B, moveset=move_set)
        # 4.1 alignment doesn't alter the underlying sequences
        @testset "4.1 alignment doesn't alter the underlying sequences" begin
            @test ungap(aligned_A) == A && ungap(aligned_B) == B
        end
        # 4.2 aligned_sequences have the same length
        @testset "4.2 aligned_sequences have the same length" begin
	        @test length(aligned_A) == length(aligned_B)
        end
        # 4.3 indels don't break codon readingframe
        @testset "4.3 indels don't break codon readingframe" begin
        # extension at edges turned off
        # reference informed alignment
        aligned_A, aligned_B = seed_chain_align(ref = A, query = B, moveset=move_set)
        # check if any clean up was needed
        @test ungap(aligned_A) == A
        @test ungap(aligned_B) == B
        end 
    end

    @testset "5. msa_codon_align" begin
        Random.seed!(42)
        m = 10
        ref = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 51)
        seqs = Vector{LongDNA{4}}(undef, m)
        for i in 1:m
            seqs[i] = mutateSequence(ref, codon_indel_avg = 1.0, frameshift_indel_avg = 1.0, sub_mutation_avg = 1.0, verbose=false)
        end
	    # 5.1 aligned_sequences have the same length
        @testset "5.1 aligned_sequences have the same length" begin
            msa = msa_codon_align(ref, seqs)
            n = length(msa[1])
            for i in 1:m
                @test n == length(msa[i+1])
            end
        end
        # 5.2 test indels not break codon readingframe after clean_up
        @testset "5.2 test indels not break codon readingframe after clean_up" begin
            msa = msa_codon_align(ref, seqs)
            cleaned_msa = clean_frameshifts(msa[1], msa[2:end])
            for i in 1:10
                @test cleaned_msa[i+1] == msa[i+1]
            end
        end
        # 5.3 type inferrence
        @testset "5.3 type inferrence" begin
            @test typeof(@inferred msa_codon_align(ref, seqs)) == Vector{LongDNA{4}}
        end
        # 5.4 scaffold_msa_from_pairwise (internal method in msa_codon_align) goes from pairwise alignments to MSA
        @testset "5.4 scaffold_msa_from_pairwise (internal method in msa_codon_align) handles codon insertions as expected" begin
            aligned_ref1 = LongDNA{4}("ATG---TTTCCCGGGTAA")
            aligned_seq1 = LongDNA{4}("ATGTTTTTTCCCGGGTAA")
            aligned_ref2 = LongDNA{4}("---ATG---TTTCCCGGGTAA---")
            aligned_seq2 = LongDNA{4}("ATGATGTTTTTTCCCGGGTAAGGG")
            aligned_ref3 = LongDNA{4}("ATGTTTCCCGGGTAA")
            aligned_seq3 = LongDNA{4}("ATG---CCCGGG---")
            msa = SeededAlignment.scaffold_msa_from_pairwise([aligned_ref1,aligned_ref2,aligned_ref3], [aligned_seq1,aligned_seq2,aligned_seq3])
            # compare with expected scaffolded msa
            @test msa[1] == LongDNA{4}("---ATG---TTTCCCGGGTAA---")
            @test msa[2] == LongDNA{4}("---ATGTTTTTTCCCGGGTAA---")
            @test msa[3] == LongDNA{4}("ATGATGTTTTTTCCCGGGTAAGGG")
            @test msa[4] == LongDNA{4}("---ATG------CCCGGG------")
        end
    end
    
	@testset "6. clean_frameshifts - pairwise" begin
        # 6.1 cleaning alignments which don't have frameshift mutations
        @testset "6.1 cleaning alignments which don't have frameshift mutations" begin
            A_no_frameshift = LongDNA{4}("---ATG---CCATTG---")
            B_no_frameshift = LongDNA{4}("TTTATGTCC---TTGTTT")
            A_cleaned, B_cleaned = clean_frameshifts(A_no_frameshift,B_no_frameshift)
            @test A_cleaned == A_no_frameshift
            @test B_cleaned == B_no_frameshift
        end
        # 6.2 clean single insertion and single deletion
        @testset "6.2 clean single insertion and single deletion" begin
            A_frameshift = LongDNA{4}("ATG-CCA")
            B_frameshift = LongDNA{4}("ATGTCC-")
            A_cleaned, B_cleaned = clean_frameshifts(A_frameshift,B_frameshift)
            @test A_cleaned == LongDNA{4}("ATGCCA")
            @test B_cleaned == LongDNA{4}("ATGCCN")
        end
        # 6.3 insertion of length 4 and deletion of length 4 which starts on codon boundary
        @testset "6.3 insertion of length 4 and deletion of length 4 which starts on codon boundary" begin
            A_frameshift = LongDNA{4}("ATG----CCATTG")
            B_frameshift = LongDNA{4}("ATGTTTT----TG")
            A_cleaned, B_cleaned = clean_frameshifts(A_frameshift,B_frameshift)
            @test A_cleaned == LongDNA{4}("ATG---CCATTG")
            @test B_cleaned == LongDNA{4}("ATGTTT---NTG")
        end
        # 6.4 insertions and deletions which cross codon boundaries and doesn't start on a boundary
        @testset "6.4 insertions and deletions which cross codon boundaries and doesn't start on a boundary" begin
            A_frameshift = LongDNA{4}("ATGA-----CCATT")
            B_frameshift = LongDNA{4}("ATGATTTTT---TG")
            A_cleaned, B_cleaned = clean_frameshifts(A_frameshift,B_frameshift)
            @test A_cleaned == LongDNA{4}("ATGACCATT")
            @test B_cleaned == LongDNA{4}("ATGANNNTG")
        end
        # 6.5 type inferrence
        @testset "6.5 type inferrence" begin
            A_frameshift = LongDNA{4}("ATGA-----CCATT")
            B_frameshift = LongDNA{4}("ATGATTTTT---TG")
            @test typeof(@inferred clean_frameshifts(A_frameshift,B_frameshift)) == Tuple{LongDNA{4},LongDNA{4}}
        end
    end

    @testset "7. clean_frameshifts - MSA" begin
        # simple test if msa version matches pairwise version
        @testset "7.1 matches pairwise version" begin
            A_frameshift = LongDNA{4}("ATGA-----CCATT")
            B_frameshift = LongDNA{4}("ATGATTTTT---TG")
            clean_msa = clean_frameshifts(A_frameshift,[B_frameshift])
            @test clean_msa[1] == LongDNA{4}("ATGACCATT")
            @test clean_msa[2] == LongDNA{4}("ATGANNNTG")
        end
        # 7.2 cleans msa example correctly
        @testset "7.2 cleans msa example correctly" begin
            aligned_seqs = Vector{LongDNA{4}}(undef, 4)
            cleaned_seqs = Vector{LongDNA{4}}(undef, 4)
            # original ref:   LongDNA{4}("ATGTTTCCCGGGTAA")
            aligned_ref =     LongDNA{4}("ATG---TTTCCCGGGT-AA")
            aligned_seqs[1] = LongDNA{4}("-TG------CCCGGGT-A-")
            aligned_seqs[2] = LongDNA{4}("ATGAAATTTCCCGGGT-AA")
            aligned_seqs[3] = LongDNA{4}("ATGAAA----CCGGGT-AA")
            aligned_seqs[4] = LongDNA{4}("ATG---TTTCCCGGGTTAA")
            # expected clean_up results
            cleaned_ref =     LongDNA{4}("ATG---TTTCCCGGGTAA")
            cleaned_seqs[1] = LongDNA{4}("NTG------CCCGGGTAN")
            cleaned_seqs[2] = LongDNA{4}("ATGAAATTTCCCGGGTAA")
            cleaned_seqs[3] = LongDNA{4}("ATGAAA---NCCGGGTAA")
            cleaned_seqs[4] = LongDNA{4}("ATG---TTTCCCGGGTAA")
            # clean and test results
            cleaned_msa = clean_frameshifts(aligned_ref, aligned_seqs)
            @test cleaned_msa[1] == cleaned_ref
            @test cleaned_msa[2] == cleaned_seqs[1]
            @test cleaned_msa[3] == cleaned_seqs[2]
            @test cleaned_msa[4] == cleaned_seqs[3]
            @test cleaned_msa[5] == cleaned_seqs[4]
        end
        # 7.3 type inferrence
        @testset "7.3 type inferrence" begin
            A_frameshift = LongDNA{4}("ATGA-----CCATT")
            B_frameshift = LongDNA{4}("ATGATTTTT---TG")
            @test typeof(@inferred clean_frameshifts(A_frameshift,[B_frameshift])) == Vector{LongDNA{4}}
        end
    end
end