# Types

This page describes the core types that are used in the Package. These types are used for repesent dna sequences, alignment movement options and scoring of alignments. For reference the core type are presented below 

## Core Types

- [`Move`](@ref) - represents an allowed alignment move during alignment (e.g. insert gap, or match nucleotide)
- [`MoveSet`](@ref) - represents a collection of `Move` objects

- [`ScoreScheme`](@ref) - represents how different aspects of an alignment is scored (e.g. substitution penalty, extension penalty)
- `LongDNA{4}` - represents the dna sequences. See [`LongDNA`](https://biojulia.dev/BioSequences.jl/stable) in the BioSequences.jl documentation.

---
In this document we outline in-depth how the `Move`, `MoveSet` and `ScoreScheme` types work. For documentation on `LongDNA{4}` we recommend looking at the documentation from [BioSequences.jl]. 

## `Move`

Represents a single allowed move during alignment (e.g., a match, insertion, or deletion). Moves can preserve codon reading frames by aligning phases.

```julia
Move(
    step::Int, 
    score::Float64,         
    horizontal_stride::Int, 
    horizontal_phase::Int,  
    vertical_stride::Int,   
    vertical_phase::Int,
    extensionAble::Bool
)
```

### Fields

- **`step::Int`**: How many alignment matrix steps this move advances. E.g. a move could correspond to a gap of length 1 or a gap of lenth 3. 
- **`score::Float64`**: The cost (or negative score) for using this move once in an alignment. If a `Move` has a high score it will be used less in the optimal alignment. 
- **`horizontal_stride::Int`**: Bases consumed from the reference sequence.
- **`horizontal_phase::Int`**: Frame offset in the reference (top sequence).
- **`vertical_stride::Int`**: Bases consumed from the query sequence.
- **`vertical_phase::Int`**: Frame offset in the query (bottom sequence).

- **`extensionAble::Bool`**: Whether this move can be extended (e.g., in gap extension scoring).

The stride and phase fields are what enable us to make some moves Codon-preserving. 

### Constructors

For convenience, there are a few constructors:

```julia
# keyword constructor
Move(; step::Int, score::Float64, stride::Int, phase::Int, extensionAble::Bool=false)
# 
Move(step, score, stride, phase, extensionAble=false)
# ignore stide and phase considerations by letting stride = 1, phase = 0
Move(step, score, extensionAble=false)
```

---

## `MoveSet`

A `MoveSet` groups together the allowable match, insertion, and deletion moves used in dynamic programming.

```julia
MoveSet(
    match_moves::Vector{Move},
    hor_moves::Vector{Move},
    vert_moves::Vector{Move}
)
```

### Fields

- **`match_moves`**: Moves that align bases from both sequences (e.g., match/mismatch).
- **`hor_moves`**: Gaps in the reference sequence (insertions).
- **`vert_moves`**: Gaps in the query sequence (deletions).

### Example

The default codon-aware move set is provided by:

```julia
std_codon_moveset()
```

An example for a custom moveset is given by
```julia
match_moves = [Move(1,.0)]
hor_moves =  [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]
vert_moves = [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]

moveset = MoveSet(match_moves=match_moves, hor_moves=hor_moves, vert_moves = vert_moves)
```

---

## `ScoreScheme`

Defines the scoring parameters used in alignment.

```julia
ScoreScheme(
    match_score::Float64,
    mismatch_score::Float64,
    extension_score::Float64,
    edge_ext_begin::Bool,
    edge_ext_end::Bool,
    kmerlength::Int
)
# default values
ScoreScheme(; match_score=0.0, mismatch_score=0.5,extension_score=0.1,edge_ext_begin=true,edge_ext_end=true,kmerlength=21)
```


### Fields

- **`match_score`**: Score for matching bases (typically 0).
- **`mismatch_score`**: Penalty for a mismatch.
- **`extension_score`**: Penalty for extending a gap.
- **`edge_ext_begin`**: Allow gap extension from the beginning of a sequence.
- **`edge_ext_end`**: Allow gap extension from the end of a sequence.
- **`kmerlength`**: Length of kmers used in alignment seeding (if applicable).

### Examples

Use the default scoring with:

```julia
scoreScheme = std_codon_scoring()
```
For custom scoring we can se the keyword constructor. Fields that are left out are kept at their default values. 
```julia
scoreScheme = scoreScheme(match_score = 0.0, mismatch_score = 0.5, extension_score = 0.3)
```
---
## LongDNA{4}

See [`LongDNA`](https://biojulia.dev/BioSequences.jl/stable) in the BioSequences.jl documentation.

---

## See Also

- [API Reference](api.md)