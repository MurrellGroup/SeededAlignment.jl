# SeededAlignment.jl

A Julia package for efficent, flexible, reading frame respecting alignments. 

## 📦 Installation

You can install this package using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/ProjectName.jl")
```

## Core Features 🔧

1. Fast flexible Codon_aligner which respects reference sequence reading frame. 
2. Fast (subquadratic) flexible pairwise_aligner, with customizable alignment moves and scoring options.
3. Clean single indel noise from Codon_Alignment

## Examples: 

Below we outline a simple usage case

```julia
using SeededAlignment

# read in dna_sequences
seq_names, dna_seqs = read_fasta("example.fasta")
# choose a reference sequence with trusted reading frame
ref_seq = dna_seqs[1]
# produce a codon alignment based on the reference sequence reading frame
alignment = codon_ref_aligner(ref_seq, dna_seqs[2:end])
# write alignment to fasta file
write_fasta("example_output.fasta", alignment, seq_names = seq_names)
```

We can customize the aligner by specifying a MoveSet and a ScoringScheme.

```


```


