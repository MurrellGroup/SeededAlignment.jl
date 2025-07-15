"""
    read_fasta(filepath::String)

Reads in a fasta file and returns a tuple of (seqnames, seqs).
"""
function read_fasta(filepath::String; top_seq_is_ref = false)
    reader = FASTX.FASTA.Reader(open(filepath, "r"))
    fasta_in = [record for record in reader]
    close(reader)
    return [String(FASTX.FASTA.identifier(rec)) for rec in fasta_in],
    [LongDNA{4}(uppercase(String(FASTX.FASTA.sequence(rec)))) for rec in fasta_in]
end

"""
    write_fasta(filepath::String, sequences::Vector{LongDNA{4}}; seq_names = nothing)

Writes a fasta file from a vector of sequences, with optional seq_names.
"""
function write_fasta(filepath::String, sequences::Vector{LongDNA{4}}; seq_names = nothing)
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
    write_fasta(filepath::String, sequences::Tuple{LongDNA{4},LongDNA{4}}; seq_names = nothing)

Writes a fasta file from a Tuple of sequences, with optional seq_names.
"""
function write_fasta(filepath::String, sequences::Tuple{LongDNA{4},LongDNA{4}}; seq_names = nothing)
    if seq_names === nothing
        seq_names = ["S$(i)" for i = 1:length(sequences)]
    end
    writer = FASTX.FASTA.Writer(open(filepath, "w"))
    for i = 1:2
        rec = FASTX.FASTA.Record(seq_names[i], sequences[i])
        write(writer, rec)
    end
    close(writer)
end