```@meta
CurrentModule = SeededAlignment
```

# API Reference

This section provides detailed documentation for all exported types and functions in `SeededAlignment.jl`.

## Core Methods

- **`seed_chain_align`** 
- **`msa_codon_align`** 
- [`clean_frameshifts`](@ref)
- [`nw_align`](@ref)

## Utilities

methods for reading and writing of sequence data in FASTA format.

- [`read_fasta`](@ref)
- [`write_fasta`](@ref)

## Core Types

- `Move`
- `Moveset`
- `ScoringScheme`
- `LongDNA{4}`

## Standard Parameters

- `std_scoring`
- `std_codon_moveset`
- `std_noisy_moveset`

## Index

A complete alphabetical listing of all documented functions.

```@index
```

## Index Docstrings

```@autodocs
Modules = [SeededAlignment]
Order = [:type, :function]
```