```@meta
CurrentModule = SeededAlignment
```

# API Reference

This section provides detailed documentation for all public types and functions in `SeededAlignment.jl`.

## Core Types

These types define the fundamental structures used in the alignment algorithms:

- `Move`
- `MoveSet`
- `ScoreScheme`
- `LongDNA{4}`

## Alignment Functions

These functions perform the actual sequence alignment operations:

- `seed_chain_align`
- `nw_align`
- `msa_codon_align`
- `clean_alignment_readingframe`

## Utilities

Helper functions for working with alignment data:

- `read_fasta`
- `write_fasta`

## Standard Parameters

Pre-configured settings for common alignment scenarios:

- `std_codon_scoring`
- `std_codon_moveset`
- `pairwise_noisy_moveset`

## Function Index

A complete alphabetical listing of all documented functions.

```@index
```

```@autodocs
Modules = [SeededAlignment]
```