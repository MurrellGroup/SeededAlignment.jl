```@meta
CurrentModule = SeededAlignment
```

# SeededAlignment

Documentation for [SeededAlignment](https://github.com/MurrellGroup/SeededAlignment.jl).

## Introduction

`SeededAlignment.jl` is a flexible high-performance Julia package for aligning DNA sequences. In particular, it can be used to create alignments that **preserves reading frame** in coding regions. The package provides tools for:

- High-performance pairwise alignments using k-mer seeds and chaining
- Producing reading frame preserved multiple sequence alignments
- Customizable alignment operations and scoring systems
- Cleaning single indel noise from pairwise codon alignment

This documentation provides a deeper look into the package’s functionality, design, and customization options.

If you're looking for installation, check out the [README on GitHub](https://github.com/MurrellGroup/SeededAlignment.jl).

---

## When to Use SeededAlignment.jl

This package is ideal for:

- Aligning **coding DNA sequences** where preserving the reading frame is critical
- Performing **fast pairwise alignments** using seed-and-extend strategies
- Creating **codon-aware multiple sequence alignments** with a trusted reference
- Customizing alignment scoring to handle mismatches, gaps, or biological constraints

---

## Core Features

- 🔗 **Seeded Pairwise Alignment** — Fast alignment using k-mer seeding and chaining  
- 🧬 **Codon-Aware MSA** — Align multiple sequences while preserving codon structure  
- ⚙️ **Custom Move Sets** — Define allowed alignment operations (e.g., frameshifts, gaps)  
- 🧠 **Flexible Scoring** — Easily adjust match/mismatch/gap penalties
- 🧹 **clean reading frame** - Remove indels which break reading frame for a pairwise alignment  
- 📄 **FASTA Support** — Simple reading and writing of sequence data  

---

## Usage Example

Here's a full example using standard codon-aware alignment tools:

```julia
using SeededAlignment
using BioSequences

# Create a reference and target sequences
ref = LongDNA{4}("ATGACGTGA")  # Reference sequence
seqs = [LongDNA{4}("ATGTCGTGA"), LongDNA{4}("ATGACGAGA")]

# Use built-in codon_alignment parameter settings
moveset = std_codon_moveset()
scoring = std_codon_scoring()

# Align sequences to the reference
alignment = msa_codon_align(ref, seqs, moveset, scoring)

# Save result to FASTA
write_fasta("alignment.fasta", alignment)
```