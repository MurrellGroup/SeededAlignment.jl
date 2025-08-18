using BioSequences
using Distributions
using Random

# TODO add an improved mutateSequence function using better noise
function mutateSequence(org_seq::LongDNA{4}; 
    codon_indel_avg::Float64 = 5.0, 
    frameshift_indel_avg::Float64 = 3.0,
    sub_mutation_avg::Float64 = 5.0    
)
    # initialize mutations
    seq = copy(org_seq)
    n = length(seq)
    dna = LongDNA{4}("ACGT")
    num_codon_indels = rand(Poisson(codon_indel_avg))
    num_frameshifts = rand(Poisson(frameshift_indel_avg))
    num_sub_mutations = rand(Poisson(sub_mutation_avg))
    # count mutation types
    cnt_codon_del = 0
    cnt_codon_ins = 0
    cnt_frm_del = 0
    cnt_frm_ins = 0
    # codon indel mutations
    for i in 1:num_codon_indels
        println(length(seq))
        indel_pos = rand(1:(n-4)÷3)
        println("indel_pos: ", 3*indel_pos)
        a = rand(Bool)
        if a 
            println("insert")
            # insertion
            # length = rand(1:4)
            len = 1
            inserted_codon = get_rand_non_stop_codons(num_codons = len)
            seq = seq[1:3*indel_pos] * inserted_codon * seq[3*indel_pos + 1 : end]
            n += 3*len
            cnt_codon_ins += len
        else
            println("deletion")
            cnt_codon_del += 1
            seq = seq[1:3*indel_pos] * seq[3*indel_pos + 4 : end]
            n -= 3
        end
        println(length(seq))
    end
    # frameshift mutations
    for i in 1:num_frameshifts
        indel_pos = rand(2:n-1)
        a = rand(Bool)
        if a
            println("frameshift insert")
            cnt_frm_ins += 1
            seq = seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[indel_pos + 1 : end]
            n += 1
        else
            # deletion
            cnt_frm_del += 1
            println("frameshift deletion")
            deleteat!(seq, indel_pos+1)
            n -= 1
        end
    end
    # substritution mutations
    for i in 1:num_sub_mutations
        println("mutation")
        mutation_pos = rand(2:n-1)
        seq = seq[1:mutation_pos-1] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[mutation_pos+1:end]
    end
    println("num codon deletions: ", cnt_codon_del)
    println("num codon insertions: ", cnt_codon_ins)
    println("num substitution mutations: ", num_sub_mutations)
    println("num frameshift deletions: ", cnt_frm_del)
    println("num frameshift insertions: ", cnt_frm_ins)
    return seq
end

function frameshift_noise_seq!(seq::LongDNA{4}; frameshift_indel_avg::Float64=3.0) # TODO Investigate why frameshift insertions look werid...
    n = length(seq)
    dna = LongDNA{4}("ACGT")
    num_frameshifts = rand(Poisson(frameshift_indel_avg))
    # insertion/deletion mutations
    for i in 1:num_frameshifts
        indel_pos = rand(1:(n-4))
        a = rand(Bool)
        if a
            # insertion
            seq = seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * seq[indel_pos + 1 : end]
            n += 1
        else
            # deletion
            deleteat!(seq, indel_pos+1)
            n -= 1
        end
    end
    return seq
end

# TODO confirm frameshifts don't become a codon indel
function frameshift_noise_seq(seq::LongDNA{4}; frameshift_indel_avg::Float64=3.0)
    n = length(seq)
    dna = LongDNA{4}("ACGT")
    num_frameshifts = rand(Poisson(frameshift_indel_avg))
    new_seq = copy(seq)
    # insertion/deletion mutations
    indel_indicies = Int64[]
    del_indicies = Int64[]
    insert_indicies = Int64[]
    for i in 1:num_frameshifts
        indel_pos = rand(2:(n-1))
        push!(indel_indicies, indel_pos)
        a = rand(Bool)
        if a
            # insertion
            push!(insert_indicies, indel_pos)
            new_seq = new_seq[1:indel_pos] * randseq(DNAAlphabet{4}(), SamplerUniform(dna"ACGT"), 1) * new_seq[indel_pos + 1 : end]
            n += 1
        else
            # deletion
            push!(del_indicies, indel_pos)
            deleteat!(new_seq, indel_pos)
            n -= 1
        end
    end
    println("num_frameshifts: ", num_frameshifts)
    println("indel_indicies: ", indel_indicies)
    println("del_indicies: ", del_indicies)
    println("insert_indicies: ", insert_indicies)
    return new_seq
end

function get_rand_non_stop_codons(; num_codons::Int64=1)
    codon_chain = LongDNA{4}("")
    non_stop_codons = get_non_stop_codons()
    for i in 1:num_codons
        append!(codon_chain, rand(non_stop_codons))
    end
    return codon_chain
end

function get_non_stop_codons()
    # bitmask represent the codons
    _DNA_BASES = ['A', 'C', 'G', 'T']
    codons = Vector{LongDNA{4}}(undef, 64)
    for bitmask in 0:63
        new_codon = LongDNA{4}("")
        # unmask the codon
        for i in 0:2
            push!(new_codon, _DNA_BASES[((bitmask >> (2 * i)) & 0x3) + 1])
        end
        codons[bitmask+1] = new_codon
    end
    non_stop_codons = filter!(x -> (BioSequences.translate(x) != LongAA('*')), codons)
    return non_stop_codons
end

function generate_random_ref(num_codons::Int64)
    ref = LongDNA{4}("ATG")
    for i in 1:(num_codons-1)
        new_codon = rand(get_non_stop_codons())
        append!(ref, new_codon)
    end
    return ref
end

#FIXME ruins the msa...
function anonymize_msa!(msa::Vector{LongDNA{4}})
    # TODO this doesn't quite work hmm...
    n = length(msa)
    for i in 1:n
        anonymize_CDS!(msa[i])
    end
end

function anonymize_CDS!(CDS::LongDNA{4})
    # anomymize original dataset
    rng = RandomDevice()
    non_stop_codons = get_non_stop_codons()
    # Note this assumes that the sequence contains no frameshifts and is CDS
    n = length(CDS) ÷ 3
    is_3gap = seq -> (seq == LongDNA{4}("---"))
    @views for i in 1:n
        cur_codon = CDS[3*(i-1)+1:3*i]
        if !is_3gap(cur_codon)
            cur_codon .= rand(rng, non_stop_codons)
        end
    end
end

#=
# TODO test if mutates CDS
function add_frameshifts!(ref_CDS::LongDNA{4}, CDS::LongDNA{4}; frameshift_indel_avg::Int64 = 3)
    @assert length(ref_CDS) == length(CDS)
    # here we allow gaps
    frameshift_indel_avg = 3
    _DNA_BASES = ['A', 'C', 'G', 'T']
    num_indels = rand(Poisson(3))
    for i in 1:num_indels
        indicies = collect(1:length(CDS))
        if rand(Bool)
            # insertion
            filter!(idx -> CDS[idx] != DNA_Gap, indicies)
            push!(indicies, length(CDS)+1)
            insert_pos = rand(indicies)
            if insert_pos == 1
                pushfirst!(CDS, rand(_DNA_BASES))
            elseif insert_pos == length(CDS)+1
                push!(CDS)
            else
                insert!(CDS, rand(_DNA_BASES),insert_pos)
            end
        else
            # deletion
            filter!(idx -> (ref_CDS[idx] != DNA_Gap && CDS[idx] != DNA_Gap), indicies)
            del_pos = rand(indicies)
            if del_pos == 1
                popfirst!(CDS)
            elseif del_pos == length(CDS)
                pop!(CDS)
            else
                CDS = CDS[1:del_pos-1] * CDS[del_pos+1:end]
            end
        end
    end
end

# TODO test if mutates CDS
function add_codon_indels!(CDS::LongDNA{4})
    @assert length(CDS) % 3 == 0
    @assert length(CDS) == length(ungap(CDS))
    # here we allow gaps
    num_codons = length(CDS) ÷ 3
    codon_indel_avg = 4.8
    _DNA_BASES = ['A', 'C', 'G', 'T']
    num_indels = rand(Poisson(codon_indel_avg))
    for i in 1:num_indels
        # TODO maybe add length of indels
        indicies = collect(1:num_codons)
        if rand(Bool) 
            # insertion
            push!(indicies, num_codons+1)
            insert_pos = rand(indicies)
            rand_codon = rand(get_non_stop_codons())
            if insert_pos == 1
                CDS = rand_codon * CDS
            elseif insert_pos == num_codons+1
                CDS = CDS * rand_codon
            else
                CDS = CDS[1:3*(idx-1)] * rand_codon * CDS[3:3*(idx-1)+1]
            end
            num_codons += 1
        else 
            # deletion
            del_pos = rand(indicies)
            if del_pos == 1
                CDS = CDS[4:end]
            elseif del_pos == num_codons
                CDS = CDS[1:length(CDS)-3]
            else
                CDS = CDS[1:3*(del_pos-1)] * CDS[3*(del_pos)+1:end]
            end
            num_codons -= 1
        end
    end
end
=#