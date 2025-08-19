import Base: show
"""

    Move(ref::Bool, step_length::Int64, score::Float64, extendable::Float64)

Represents a gap move during alignment. The `Move` instance represents either an insertion or deletion depending
on how it is passed to the `Moveset` instance - collection of `Move` instances used in alignment methods. 

For example, if Move.step_length = 3 then that represents a gap of length 3 (either insertion or deletion). 

# Extended Help

# Fields
-`ref::Bool`: Whether move respects the coding reading frame
-`step_length::Int64`: Length of gap in alignment (must be  1,2 or 3)
-`score::Float64`: penalty for using the move - cost of adding the gaps to the alignment (must be < 0)
-`extendable::Bool`: Whether the move can be affinely extended by another move. 

# Constructors
-`Move(; ref=false::Bool, step_length::Int64, score::Float64, extendable=false)`:
keyword constructor for more easily constructing Moves. Has some default values but requires at least `step_length` and `score` to be provided. 

-`RefMove(; score::Float64)`: 
Constructor for a `Move` that respect coding reading frame. Only requires `score` argument. 
Returns: `Move(ref=true, step_length=3, score=score, extendable=true)`

-`FrameshiftMove(; step_length::Int64, score::Float64, extendable::Bool=false)`: 
Constructor for a `Move` that cause frameshifts and breaks reading frame symmetry. Requires arguments `step_length`, `score`, and `extendable`.
Returns: `Move(ref=false, step_length=step_length, score=score, extendable=extendable)`

# Examples

1. codon moveset with no frameshifts

```julia
# represents codon insertion or codon deletion
codon_indel = RefMove(score=1.0)
#= passing to moveset solidifies what the allowed alignment operations are,
namely single codon insertions and deletions. 
=#
Moveset(ref_insertions = (codon_indel,), ref_deletions = (codon_indel,))
```

2. codon moveset with frameshift moves allowed

```julia
# represents codon insertion or codon deletion
codon_indel = RefMove(score=-1.0)
# represents frameshift causing insertion or deletion
frm_indel = FrameshiftMove(score=-1.5, step_length=1, extendable=true)
#= passing to moveset solidifies what the allowed alignment operations are,
namely single codon insertions and deletions and single nucleotide indels that are extendable. 
=#
Moveset(ref_insertions = (codon_indel,frm_indel), ref_deletions = (codon_indel, frm_indel))
```

"""
struct Move
    ref::Bool
    step_length::Int64
    score::Float64
    extendable::Bool
    # force construction invariants
    function Move(ref::Bool,  step_length::Int64, score::Float64, extendable::Bool)
        (step_length in (1, 2, 3)) || throw(ArgumentError("step_length must be 1, 2, or 3"))
        (score < 0) || throw(ArgumentError("score must be negative"))
        (!ref || (step_length == 3 && extendable)) || throw(ArgumentError("Invalid ref move:\nWhen ref=true it is required that step_length=3 and extendable=true."))
        new(ref, step_length, score, extendable)
    end
end

# useful constructors
function Move(; ref=false::Bool, step_length::Int64, score::Float64, extendable=false) 
    Move(ref, step_length, score, extendable)
end
# assumes ref true
function RefMove(; score::Float64) 
    Move(ref=true, step_length=3, score=score, extendable=true)
end
# assumes ref false
function FrameshiftMove(; step_length::Int64, score::Float64, extendable::Bool=false)
    Move(ref=false, step_length=step_length, score=score, extendable=extendable)
end

"""
    Moveset{X,Y}(; ref_insertions::NTuple{X,Move}, ref_deletions::NTuple{Y,Move})

Represents a collection of `Move` instances that are either insertions or deletions

# Extended Help

# Fields
-`vert_moves::NTuple{X,Move}`: insertions relative to reference - gap operation in top/first provided sequence 
-`hor_moves::NTuple{Y,Move}`: deletions relative to reference - gap operations in left/second provided sequence 

# Examples

1. codon moveset with no frameshifts

```julia
# represents codon insertion or codon deletion
codon_indel = RefMove(score=1.0)
#= passing to moveset solidifies what the allowed alignment operations are,
namely single codon insertions and deletions. 
=#
ms = Moveset(ref_insertions = (codon_indel,), ref_deletions = (codon_indel,))
```

2. codon moveset with frameshift moves allowed

```julia
# represents codon insertion or codon deletion
codon_indel = RefMove(score=-1.0)
# represents frameshift causing insertion or deletion
frm_indel = FrameshiftMove(score=-1.5, step_length=1, extendable=true)
#= passing to moveset solidifies what the allowed alignment operations are,
namely single codon insertions and deletions and single nucleotide indels that are extendable. 
=#
ms = Moveset(ref_insertions = (codon_indel,frm_indel), ref_deletions = (codon_indel, frm_indel))
```

"""
struct Moveset{X,Y}
    # gap operation in top/first provided sequence
    vert_moves::NTuple{X,Move}
    # gap operations in left/second provided sequence
    hor_moves::NTuple{Y,Move}
end
# reference informed Moveset constructor
function Moveset(; ref_insertions::NTuple{X,Move}, ref_deletions::NTuple{Y,Move}) where {X,Y}
    Moveset(ref_insertions, ref_deletions)
end
# non-reference informed Moveset constructor - alignment operations are symmetric
Moveset(gap_moves) = Moveset(gap_moves, gap_moves)

function show(io::IO, ms::Moveset)
	println(io, "Moveset(\n","  vert_moves=$(ms.vert_moves)\n","  hor_moves=$(ms.hor_moves)\n)")
end

function contains_ref_move(moveset::Moveset)
    ref_move_found = any(move.ref for move in moveset.vert_moves) || any(move.ref for move in moveset.hor_moves)
    return ref_move_found
end

# default movesets
"""
    STD_CODON_MOVESET

Constants that represents the default codon moveset with frameshift moves allowed

# default parameter values

```julia
const STD_CODON_MOVESET = Moveset(
    (
        Move(ref=false, step_length=1, score=-2.0, extendable=true),
        Move(ref=true,  step_length=3, score=-1.0, extendable=true)
    )
)
```
"""
const STD_CODON_MOVESET = Moveset(
    (
        Move(ref=false, step_length=1, score=-2.0, extendable=true),
        Move(ref=true,  step_length=3, score=-1.0, extendable=true)
    )
)
# deletions same as insertions
const STD_NOISY_MOVESET = Moveset(
    (
        Move(ref=false, step_length=1, score=-2.0, extendable=false),
        Move(ref=false, step_length=3, score=-1.0, extendable=true)
    )
)