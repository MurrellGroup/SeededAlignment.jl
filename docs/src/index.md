## SeededAlignment.jl
```@meta
CurrentModule = SeededAlignment
```

**`SeededAlignment`** provides an interface for **codon alignment** methods that use **seeding** and robustly handle **frameshift errors**. 
 
This documentation aims to more comprehensively than [README on GitHub](https://github.com/MurrellGroup/SeededAlignment.jl) explain the purpose, interface and usage of the package.

# Why use SeededAlignment?
Reading frame perserved codon alignment are essential for many modern applications such as selection detection, phylogenetic inference, and evolutionary modeling. A common problem with preparing codon alignments for such analyses is that the contig reads might contain frameshift errors due to e.g. sequencing or annotation errors. The present of such errors breaks readingframe symmetry and significantly hurts downstream analyses.

**`SeededAlignment`** removes this obstacle by detecting and removing frameshift errors during alignment; thus enabling a smoother pipeline. At the same time it also tries to provide be sufficently flexible and performant alignment methods for most biological applications. 

**Note:** a core assumption that is made by the package is that one of the input sequences has an intact trusted reading frame. This sequence is refered to as the *reference* and determines the reading frame coordinates for the other sequences in the alignment, which allows us to detect and correct frameshift errors present in the remaining sequences. 

#  Getting Started
The principal type used in the package is `LongDNA{4}` from [BioSequences.jl](https://biojulia.dev/BioSequences.jl/stable/). In essence, it implements a compact vector representation of a DNA sequence. All alignment methods only accept DNA sequences represented `LongDNA{4}`. To convert raw data into this form we can either read from a fasta file with `read_fasta`
```julia
# here seqs is of type Vector{LongDNA{4}} and seq_names is of type Vector{String}
seq_names, seqs = read_fasta("example_file.fasta")
```
or by casting a string to `LongDNA{4}`
```julia
dna_string = "ATGCTCTAA"
# this casts the string above to LongDNA{4}
LongDNA_rep = LongDNA{4}("ATGCTCTAA")
```

---

## Pairwise Alignment
To produce a pairwise codon alignment we can use the `seed_chain_align` method. There are two methods named `seed_chain_align`. The **1st** of which is for codon alignment and takes a reference sequence (sequence with trusted reading frame) and a query sequence (sequence which may or may not contain frameshift errors). An example of it is provided below.  

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("codon_example.fasta")
# select the reference sequence
ref_seq = dna_seqs[1]
# select the query sequence to be aligned against the reference
query_seq = dna_seqs[2]
# compute codon alignment with ref_seq as refernce
codon_alignment = seed_chain_align(ref=ref_seq, query=query_seq)
```

**Warning:** it is important to note that if `ref_seq` (i.e. the reference sequence) does in fact contain frameshift errors then
real atrifacts of the intended alignment might end up being corrupted. 

The **2nd** `seed_chain_align` method is for purely semantic nucleotide alignments, which we choose to refer to as De-Novo alignment. De-novo alignment treats the sequences as equals and are therefore supplied as positional arguments. An example is provided below.  

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("DeNovo_example.fasta")
# choose to sequences from the dataset
seq1 = dna_seqs[1]
seq2 = dna_seqs[2]
# align them without regard for codon boundaries
DeNovo_nucleotide_alignment = seed_chain_align(seq1, seq2)
# NOTE: Do not confuse with codon_alignment! 
```
**Caution:** Be careful not to mix these two `seed_chain_align` method up for disasterous consequences. 

**Note:** that the `nw_align` has exactly the same interface as `seed_chain_align`. In fact, in all of the previous examples `seed_chain_align` can be replaced with `nw_align` without issue.

---

## Cleaning Frameshifts
Sometimes we might prefer to clean a codon alignment that contains frameshifts as opposed to recomputing an entire alignment which could be costly. The `clean_frameshifts` method provides the means to do this. It can clean pairwise (or multiple) sequence alignment of frameshift errors and make the alignment reading frame perserved. To clean an alignment we must once again select a reference sequence with preserved reading frame. To clean a pairwise alignment we do the following

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("unclean_pairwise_alignment.fasta")
# select the aligned reference sequence
aligned_ref = dna_seqs[1]
# select the other aligned sequence
aligned_seq = dna_seqs[2]
# clean the pairwise alignment of frameshifts relateive to reference frame
cleaned_alignment = clean_frameshifts(aligned_ref, aligned_seq, verbose=true)
```
where `verbose=true` enables printout which describe what the clean up removed and added to the non-reference sequence. This can
be useful if one wants to quickly inspect what changes were made to the original sequences.

For multiple sequence alignments we can do something similar, namely
```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("unclean_msa.fasta")
# select the aligned reference sequence
aligned_ref = dna_seqs[1]
# select all the other aligned sequence
aligned_seqs = dna_seqs[2:end]
# clean the pairwise alignment of frameshifts relateive to reference frame
cleaned_msa = clean_frameshifts(aligned_ref, aligned_seqs, verbose=true)
```
where once again `verbose=true` enables descriptive printouts for the changes made. 

This latter method is particularly useful in a context where we already have computed a large multiple sequence alignment 
that contains some frameshift which we want to clean up without recomputing an expensive alignment.

---

## Codon Multiple Sequence Alignment
To produce a codon multiple sequence alignment one can use the `msa_codon_align`. First one must choose a reference sequence with preserved reading frame. The remaining sequences are then pairwise aligned against the reference - these alignments are then scaffolded into a multiple sequence alignment using the reference's frame. We show how to use the method below. 

```julia
using SeededAlignment

seq_names, dna_seqs = read_fasta("msa_example.fasta")
# choose a reference sequence from the dataset
ref_seq = dna_seqs[1]
# collect the remaining sequences to be aligned
non_ref_sequences = dna_seqs[2]
# compute codon multiple sequence alignment
codon_msa = msa_codon_alignment(ref_seq, non_ref_sequences) 
```

It is important to note that the resulting alignment will be biased by the choosen reference. This is because 
there are no direct intra sequence comparisons for the non-reference sequences.

If one believes that the reference sequence biases the codon alignment too heavily one could run this as a denoising step for detecting 
and removing frameshift errors, whereafter one could run the ungapped alignment through a different multiple sequence alignment program.
This improves the likelihood of the resulting alignment being reading frame preserved as frameshifts have been removed.

---

## API

[API Reference](api.md)
Full list of exported types, functions and constants.