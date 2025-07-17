import Base: show
"""
    Move

Represents a single allowed alignment step between two sequences in the alignment matrix.

# Fields
- `step::Int`: How many aligned units (nucleotides) this move advances by.
- `score::Float64`: Cost or score associated with performing this move.
  
### Reference to the top sequence (horizontal axis)
- `horizontal_stride::Int`: How far this move advances along the top (reference) sequence.
- `horizontal_phase::Int`: Frame offset in the top sequence; helps track codon boundaries.

### Reference to the bottom sequence (vertical axis)
- `vertical_stride::Int`: How far this move advances along the bottom (query) sequence.
- `vertical_phase::Int`: Frame offset in the bottom sequence.

- `extensionAble::Bool`: Whether this move can be extended (e.g., for affine gap penalties).

# Description
Moves define basic operations used in dynamic programming alignment: matches, mismatches, and gaps.
Moves can preserve codon reading frames by keeping strides in multiples of 3 and matching phase positions.

Two helper constructors are provided:
- `Move(step, score, stride, phase, extensionAble)` assumes reference (vertical) stride and phase
- `Move(step, score, extensionAble)` for simple moves with unit stride and zero phase.
"""
struct Move # TODO remove vertical_stride and phase for consistency
    step::Int
    score::Float64

    # refers to readingFrame of top sequence 
    horizontal_stride::Int
    horizontal_phase::Int

    # refers to readingFrame of bottom sequence
    vertical_stride::Int
    vertical_phase::Int

    # allows moves to be extended from 
    extensionAble::Bool
end

#TODO fix the useful constructor.

Move(step::Int64, score::Float64, stride::Int64, phase::Int64, extensionAble::Bool=false) =
    Move(step, score, 1,0, stride, phase, extensionAble)

#function Move(step::Int64, score::Float64, h_stride::Int64, h_phase::Int64, v_stride::Int64, v_phase::Int64, extensionAble::Bool=false)
#   return Move(step, score, h_stride, h_phase, v_stride, v_phase, extensionAble)
#end

function Move(; step::Int64, score::Float64, stride::Int64, phase::Int64, extensionAble::Bool=false)
    return Move(step, score, stride, phase, extensionAble)
end

Move(step::Int64, score::Float64, extensionAble::Bool=false) = Move(step, score, 1, 0, extensionAble)

function show(io::IO, m::Move)
	print(io, "Move(",
			  "step=$(m.step), ",
			  "score=$(m.score), ",
              "h_stride=$(m.horizontal_stride), ",
			  "h_phase=$(m.horizontal_phase), ",
			  "v_stride=$(m.vertical_stride), ",
			  "v_phase=$(m.vertical_phase), ",
			  "extensionAble=$(m.extensionAble))")
end

"""
    Moveset

Defines the full set of allowed alignment moves used during pairwise or multiple sequence alignment.

# Fields
- `vert_moves::Vector{Move}`: Moves that introduce gaps in the **query (vertical)** sequence. 
- `hor_moves::Vector{Move}`: Moves that introduce gaps in the **reference (horizontal)** sequence.

# Description
A `Moveset` groups the allowable `Move`s into categories used during dynamic programming alignment. Each type controls how the algorithm can transition between states, including nucleotide-level match and gap moves with customizable reading frame behavior.

Used by alignment algorithms such as `seed_chain_align` and  `msa_codon_align` to control the scoring and allowed operations during alignment.
"""

# stack allocated
struct Moveset{X,Y}
    vert_moves::NTuple{X,Move}
    hor_moves::NTuple{Y,Move}
end


Moveset(; vert_moves, hor_moves) = Moveset(vert_moves, hor_moves)

function show(io::IO, ms::Moveset)
	print(io, "Moveset(\n",
            "  vert_moves=$(ms.vert_moves)\n)",
            "  hor_moves=$(ms.hor_moves) \n")
end

"""
    std_codon_moveset()

Return a standard codon-aware `Moveset` for sequence alignment.

# Returns
- `Moveset`: Contains codon-aware match and gap moves.

# Moves
- Match: `Move(1, 0.0)`
- Horizontal (gaps in reference sequence):
  - `Move(1, 2.0, 1, 0, 1, 0, false)`
  - `Move(3, 2.0, 1, 0, 3, 0, true)`
- Vertical (gaps in query sequence):
  - `Move(1, 2.0, 1, 0, 1, 0, false)`
  - `Move(3, 2.0, 1, 0, 3, 0, true)`

Use this as a default `Moveset` for codon-preserving alignments.
"""
function std_codon_moveset()
    vert_moves = (Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true))
    hor_moves =  (Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true))
    return Moveset(vert_moves,hor_moves)
end
"""
    pairwise_noisy_moveset()

"""
function std_noisy_moveset()
    vert_moves = (Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,1,0, true))
    hor_moves =  (Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,1,0, true))
    return Moveset(vert_moves,hor_moves)
end