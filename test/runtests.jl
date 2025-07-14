using Test

# test dependencies
using SeededAlignment
using BioSequences
using Random

# TODO compare fast_translate to BioSequences.translate

@testset "SeededAlignment.jl" begin
    @testset "1. noisy alignment nw_affine" begin
        Random.seed!(42)
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 2001)
        B = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1980)
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
            score_params_on = ScoringScheme(match_score=0.0, mismatch_score=0.3, extension_score=0.4, edge_ext_begin=true,edge_ext_end=true)
            score_params_off = ScoringScheme(match_score=0.0, mismatch_score=0.3, extension_score=0.4, edge_ext_begin=false,edge_ext_end=false)
            move_set = Moveset(
                match_moves = [Move(1, 0.0, 1, 0, 1, 0, false)],
                hor_moves = [Move(3, 30, 1, 0, 1, 0, true)],
                vert_moves = [Move(3, 30, 1, 0, 1, 0, true)]
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

    # TODO finish this test case important
    #@testset "2. ref alignment nw_affine" begin
    #     move_set = Moveset(
    #            match_moves = [Move(1, 0.0, 1, 0, 1, 0, false)],
    #            hor_moves = [Move(3, 30, 1, 0, 1, 0, true)],
    #            vert_moves = [Move(3, 30, 1, 0, 1, 0, true)]
    #        )
    #    Random.seed!(42)
    #    A = LongDNA{4}()
    #    B = LongDNA{4}()
    #    # reference informed alignment
    #    aligned_A, aligned_B = nw_align(ref = A, query = B, clean_up_flag=true)
    #    # 2.1 indels don't break codon readingframe
    #    @test aligned_A == cleaned_A
    #    @test aligned_B == cleaned_B
    #end

    @testset "3. noisy alignment seed_chain_align" begin 
        Random.seed!(42)
        A = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 2001)
        B = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1980)
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

    # TODO finish this test case 
    #@testset "4. ref alignment seed_chain_align" begin
    #    Random.seed!(42)
    #    A = LongDNA{4}()
    #    B = LongDNA{4}()
    #    # reference informed alignment
    #    aligned_A, aligned_B = seed_chain_align(ref = A, query = B)
    #end

    @testset "5. msa_codon_align" begin
        Random.seed!(42)
        ref = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 51)
        seqs = Vector{LongDNA{4}}(undef, 20)
        for i in 1:20
            seqs[i] = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 33+i)
        end

	    # 5.1 aligned_sequences have the same length
        @testset "5.1 aligned_sequences have the same length" begin
            msa = msa_codon_align(ref, seqs)
            n = length(msa[1])
            for i in 1:20
                @test n == length(msa[i+1])
            end
        end
        # 5.2 test indels not break codon readingframe after clean_up
        @testset "5.2 test indels not break codon readingframe after clean_up" begin
            msa = msa_codon_align(ref, seqs)
            aligned_ref = msa[1]
            for i in 1:20
                cleaned_ref, cleaned_seq = clean_frameshifts(aligned_ref, msa[i+1])
                @test cleaned_seq == msa[i+1]
            end
        end
        # 5.3 type inferrence
        @testset "5.3 type inferrence" begin
            @test typeof(@inferred msa_codon_align(ref, seqs)) == Vector{LongDNA{4}}
        end

    end
    
	@testset "6. clean_frameshifts" begin
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
end