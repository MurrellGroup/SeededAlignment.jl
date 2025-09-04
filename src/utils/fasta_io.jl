"""
    read_fasta(filepath::String)

Reads in a fasta file and returns a tuple of (seqnames, seqs). 
"""
    function read_fasta(filepath::String)
        reader = FASTX.FASTA.Reader(open(filepath, "r"))
        fasta_in = [record for record in reader]
        close(reader)
        seq_names = [String(FASTX.FASTA.identifier(rec)) for rec in fasta_in]
        try
            dna_seqs = [LongDNA{4}(uppercase(String(FASTX.FASTA.sequence(rec)))) for rec in fasta_in]
            return seq_names, dna_seqs
        catch
            try 
                aa_seqs = [LongAA(uppercase(String(FASTX.FASTA.sequence(rec)))) for rec in fasta_in]
                return seq_names, aa_seqs
            catch
                string_seqs = [String(FASTX.FASTA.sequence(rec)) for rec in fasta_in]
                return seq_names, string_seqs
            end
        end
    end

"""
    write_fasta(filepath::String, sequences::Union{Vector{LongDNA{4}}, Vector{LongAA}, Vector{String}}; seq_names = nothing)

Writes a fasta file from a vector of sequences, with optional seq_names.
"""
function write_fasta(filepath::String, sequences::Union{Vector{LongDNA{4}}, Vector{LongAA}, Vector{String}}; seq_names = nothing)
    if seq_names === nothing
        seq_names = ["S$(i)" for i = 1:length(sequences)]
    end
    writer = FASTX.FASTA.Writer(open(filepath, "w"))
    for i = 1:length(seq_names)
        rec = FASTX.FASTA.Record(seq_names[i], sequences[i])
        write(writer, rec)
    end
    close(writer)
end

"""
    write_fasta(filepath::String, sequences::Union{NTuple{N,LongDNA{4}}, NTuple{N,LongAA}, NTuple{N,String}}; seq_names = nothing)

Writes a fasta file from a Tuple of sequences, with optional seq_names.
"""
function write_fasta(filepath::String, sequences::Union{NTuple{N,LongDNA{4}}, NTuple{N,LongAA}, NTuple{N,String}}; seq_names = nothing) where {N}
    if seq_names === nothing
        seq_names = ["S$(i)" for i = 1:length(sequences)]
    end
    writer = FASTX.FASTA.Writer(open(filepath, "w"))
    for i = 1:length(sequences)
        rec = FASTX.FASTA.Record(seq_names[i], sequences[i])
        write(writer, rec)
    end
    close(writer)
end