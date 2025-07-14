# SeededAlignment.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MurrellGroup.github.io/SeededAlignment.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/SeededAlignment.jl/actions/workflows/CI.yml?query=branch%3Amain)

SeededAlignment.jl aims to provide a user-friendly interface for pairwise sequence alignment that supports user-defined alignment operations, flexible scoring and is optimized for performance. In particular, it can be used to produce frameshift-free codon-alignments (relative to known refence) and clean previously constructed alignments of frameshift mutations, provided an aligned reference sequence.

## Core Methods

1. **`seed_chain_align`** - fast (in the sense of subquadratic TC) pairwise alignment method which uses heurstically guided seeding and subsequently chains them together.  
2. **`msa_codon_align`** - scaffolds a frameshift-free multiple sequence alignment from pairwise alignments relative to a provided reference with intact reading frame.  
3. **`clean_frameshifts`** - clean pairwise alignments of frameshift mutations given that one of the sequences has an intact reading frame.  
4. **`nw_align`** - classic Needleman-Wunsch algorithm for pairwise alignment.

**FASTA Support** – simple reading and writing of sequence data in FASTA format.

--- 

## Some Use Cases in Production Pipeline

1. **Sequence Clustering**  
   Leveraging `seed_chain_align` to define custom distance or similarity metrics for clustering sequences. This can be used for simple exploratory analysis or if one wants to perform MSA using clustering information as a guide.

2. **Quality Control – Sequencing and Annotation Errors**  
   Use `seed_chain_align` as part of a quality control pipeline to detect sequencing or annotation issues such as outliers, chimeras, or contamination by aligning sequences against each other.
 
3. **Clean problematic alignments that contain frameshift mutations, (provided one of the sequences has intact reading frame)**
   Often, we have a sequence alignment that contains frameshift mutations due to sequencing or annotation errors. These errors can be "corrected" automatically by using the `clean_frameshifts` method to produce a frameshift-free codon pairwise alignments (support for MSA with aligned_reference coming soon).

4. **Produce frameshift-free visual multiple sequence alignment (based on given reference sequence).**  
   Many downstream applications, such as selection analysis and phylogenetic inference, require frameshift-free multiple sequence alignments as part of their model assumptions. Provided a reference sequence with intact reading frame, one can use `msa_codon_align` to get a visual frameshift-free multiple sequence alignment based on the pairwise alignments relative to the reference. 

---

## Installation

You can install this package using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/MurrellGroup/SeededAlignment.jl")
```

---

## Examples:

Here we provide quick showcase of the following methods:

- `seed_chain_align` - with and without refernce
- `clean_frameshifts` - given pairwise alignment with reference. 
- `msa_codon_align` - produces a visual multiple sequence alignment pairwise relative to a reference

### seed_chain_align 

seed_chain_align is provided with two wrapper methods for two different uses. The first of which looks like
```julia
seed_chain_align(seq1, seq2)
```
and makes no assumptions about the two sequences. On the other hand if we have a reference sequence with intact readingframe we have to specify this to the method. This is done by calling the other wrapper which looks like
```julia
seed_chain_align(ref=ref_seq, query=non_ref_seq)
```
the main difference between the results of two methods is that the latter prefers codon indels relative to reference readingframe over indels which shift the readingframe.

**Note:** that the nw_align wrappers are the same as seed_chain_align. Hence in all of the following examples seed_chain_align can be replaced with nw_align without any issues. 

If we want to align two sequnces without assuming intact readingframe we simply supply them as positional arguments. 

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("example.fasta")
# choose two sequences to align
seq1 = dna_seqs[1]
seq2 = dna_seqs[2]
# align heurisitically where no sequences is treated as reference
seq1_aligned, seq2_aligned = seed_chain_align(seq1,seq2)
```
On the other hand if we have a reference sequence with intact readingframe we have to specify this to the method

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("example.fasta")
# choose two sequences to align
ref_seq = dna_seqs[1]
non_ref_seq = dna_seqs[2]
# align heurisitically where ref_seq is treated as reference
seq1_aligned, seq2_aligned = seed_chain_align(ref=ref_seq, query=non_ref_seq)
```

### clean_frameshifts 

If we want to make sure that we get a valid codon alignment (on the nucleotide level) with no frameshift mutations - then we can remove the frameshift mutations via calling the `clean_frameshifts` method on the final result. This can be done in two ways: either via supplying a boolean argument to the alignment method or by manually calling `clean_frameshifts` yourself. This is shown below. 

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("example.fasta")
# choose two sequences to align
ref_seq = dna_seqs[1]
non_ref_seq = dna_seqs[2]
# align heurisitically where ref_seq is treated as reference and cleans up frameshift mutations from final alignment
seq1_aligned, seq2_aligned = seed_chain_align(ref=ref_seq, query=non_ref_seq, clean_up_enabled=true)
```
or this can be done manually
```julia
# align heurisitically where ref_seq is treated as reference
seq1_aligned, seq2_aligned = seed_chain_align(ref=ref_seq, query=non_ref_seq)
# cleans up frameshift mutations
cleaned_seq1, cleaned_seq2 = clean_frameshifts(seq1_aligned, seq2_aligned)
```
**Note:** that one can supply a verbose flag for visibility of what edits were made during clean up:
```julia
seq1_aligned, seq2_aligned = seed_chain_align(ref=ref_seq, query=non_ref_seq, clean_up_enabled=true, verbose=true)
cleaned_seq1, cleaned_seq2 = clean_frameshifts(seq1_aligned, seq2_aligned, verbose=true)
```

### msa_codon_align

Lastly, we show how `msa_codon_align` can be used to visual frameshift-free multiple sequence alignment based on pairwise alignments relative to a reference. Note that some edits to the sequences might be made since we clean framshift mutations in each pairwise alignment. 

**NOTE:** A refernce sequence with intact readingframe is a required argument for this method to work.

```julia
using SeededAlignment

# read in dna_sequences
seq_names, dna_seqs = read_fasta("example.fasta")
# choose a reference sequence with intact readingframe
ref_seq = dna_seqs[1]
# extract the query sequences
query_seqs = dna_seqs[2:end]
# produce a codon alignment based on the reference sequence reading frame
alignment = msa_codon_align(ref_seq, query_seqs)
# write alignment to fasta file
write_fasta("example_output.fasta", alignment, seq_names = seq_names)
```