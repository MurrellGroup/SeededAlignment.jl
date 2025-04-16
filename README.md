# SeededAlignment.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/SeededAlignment.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for efficent, flexible, reading frame respecting alignments. 

## 📦 Installation

You can install this package using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/MurrellGroup/SeededAlignment.jl")
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
moveset = std_codon_moveset()
scoreScheme = std_codon_scoring()
# produce a codon alignment based on the reference sequence reading frame
alignment = msa_codon_align(ref_seq,dna_seqs[2:end], moveset, scoreScheme)
# write alignment to fasta file
write_fasta("example_output.fasta", alignment, seq_names = seq_names)
```

We can specifying custom a MoveSet and ScoringScheme to guide the alignment.

```Julia
using SeededAlignment

# read in dna_sequences
seq_names, dna_seqs = read_fasta("example.fasta")
# choose a reference sequence with trusted reading frame
ref_seq = dna_seqs[1]

# customizing moveset
match_moves = [Move(1,.0)]
hor_moves =  [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]
vert_moves = [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]
moveset = MoveSet(match_moves,hor_moves,vert_moves)

# customizing scoreScheme
score_params = ScoreScheme(match_score = -0.5, mismatch_score = 1.3 extension_score = 0.7, kmerlength = 24)
# NOTE: left out fields are kept at default values. Check documentation of ScoreScheme to see default values. 

# produce a codon alignment based on the reference sequence reading frame
alignment = msa_codon_align(ref_seq, dna_seqs[2:end], moveset, score_params)
# write alignment to fasta file
write_fasta("example_output.fasta", alignment, seq_names = seq_names)
```