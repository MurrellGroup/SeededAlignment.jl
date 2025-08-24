```@meta
CurrentModule = SeededAlignment
```

# API Reference

This section provides detailed documentation for all exported types, functions and constants in `SeededAlignment.jl`.

## Core Methods

- [`seed_chain_align`](@ref)
- [`msa_codon_align`](@ref) 
- [`clean_frameshifts`](@ref)
- [`nw_align`](@ref)

## Utilities

methods for reading and writing of sequence data in FASTA format.

- [`read_fasta`](@ref)
- [`write_fasta`](@ref)

## Core Types

- [`Move`](@ref)
- [`Moveset`](@ref)
- [`ScoringScheme`](@ref)
- [`LongDNA{4}`] - Array specialized for DNA from BioSequences.jl

## Constants

- [`STD_SCORING`](@ref)
- [`STD_CODON_MOVESET`](@ref)
- [`STD_NOISY_MOVESET`](@ref)

## Index

A complete alphabetical listing of all documented functions.

```@index
```

## Index Docstrings

```@autodocs
Modules = [SeededAlignment]
Order = [:type, :function, :constant]
```