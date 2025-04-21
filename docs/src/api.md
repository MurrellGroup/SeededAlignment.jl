```@meta
CurrentModule = SeededAlignment
```

# API Reference

This section provides detailed documentation for all public types and functions in `SeededAlignment.jl`.

## Core Types

These types define the fundamental structures used in the alignment algorithms:

- [`Move`](@ref)
- [`MoveSet`](@ref)
- [`ScoreScheme`](@ref)
- `LongDNA{4}`

## Alignment Functions

These functions perform the actual sequence alignment operations:

- [`seed_chain_align`](@ref)
- [`msa_codon_align`](@ref)
- [`clean_alignment_readingframe`](@ref)
- [`nw_align`](@ref)

## Utilities

Helper functions for working with alignment data:

- [`read_fasta`](@ref)
- [`write_fasta`](@ref)

## Standard Parameters

Pre-configured settings for common alignment scenarios:

- [`std_codon_scoring`](@ref)
- [`std_codon_moveset`](@ref)
- `pairwise_noisy_moveset`

## Function Index

A complete alphabetical listing of all documented functions.

```@index
```

```@autodocs
Modules = [SeededAlignment]
Order = [:type, :function]
```