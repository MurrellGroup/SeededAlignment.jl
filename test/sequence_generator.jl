function generate_seq_pair(seq_len::Int64, long_indel_frequency::Float64, triple_indel_frequency::Float64, substitution_frequency::Float64, frameshift_frequency::Float64, long_length::Int64)
    seq = randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), seq_len)
    dna = LongDNA{4}("ACGT")
   
    # mutations

    mutated_seq = LongDNA{4}()

    for i in seq
        push!(mutated_seq, i)
    end

    n = seq_len

    triple_indels = rand(Poisson(n * triple_indel_frequency))
    frameshifts = rand(Poisson(n * frameshift_frequency))
    long_indels = rand(Poisson(n * long_indel_frequency))
    substitutions = rand(Poisson(n * substitution_frequency))

    for i in 1:substitutions
        mutated_seq[rand(1:n)] = dna[rand(findall(x->x != mutated_seq[i], dna))]
    end

    for i in triple_indels
        if rand(1:2) == 1

            # Insertion
            insertion_pos = rand(1:n)
            if insertion_pos < n
                mutated_seq = mutated_seq[1:insertion_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 3) * mutated_seq[insertion_pos + 1:end]
            else
                mutated_seq = mutated_seq[1:insertion_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 3)
            end
        else
            # deletion

            deletion_pos = rand(1:n - 2)
            for i in 1:3
                deleteat!(mutated_seq, deletion_pos)
            end
        end
        n = length(mutated_seq)
    end

    for i in 1:frameshifts
        indel_pos = rand(1:n)
        if rand(1:2) == 1
            # insertion
            mutated_seq = mutated_seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * mutated_seq[indel_pos + 1 : end]
        else
            #deletion
            deleteat!(mutated_seq, indel_pos)
        end
        n = length(mutated_seq)
    end

    for i in 1:long_indels
        len = rand(Poisson(long_length/3)) * 3

        if rand(1:2) == 1
            # insertion
            insertion_pos = rand(1:n)
            if insertion_pos < n
                mutated_seq = mutated_seq[1:insertion_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), len) * mutated_seq[insertion_pos + 1:end]
            else 
                mutated_seq = mutated_seq[1:insertion_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), len)
            end
        else
            # delete
            deletion_pos = rand(1:n + 1 - len)
            for i in 1:len
                deleteat!(mutated_seq, deletion_pos)
            end
        end
        n = length(mutated_seq)
    end

    return seq, mutated_seq
end

#Example: generate_seq_pair(50, 0.2, 0.1, 0.01, 0.01, 4)