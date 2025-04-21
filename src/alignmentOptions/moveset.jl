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

Move(; step::Int64, score::Float64, stride::Int64, phase::Int64, extensionAble::Bool=false) = 
    Move(step::Int64, score::Float64, stride::Int64, phase::Int64, extensionAble::Bool)

Move(step::Int64, score::Float64, stride::Int64, phase::Int64, extensionAble::Bool=false) =
    Move(step, score, stride, phase, 1, 0, extensionAble)

Move(step::Int64, score::Float64, extensionAble::Bool=false) = Move(step, score, 1, 0, extensionAble)

function show(io::IO, m::Move)
	print(io, "Move(",
			  "step=$(m.step), ",
			  "score=$(m.score), ",
			  "v_stride=$(m.vertical_stride), ",
			  "v_phase=$(m.vertical_phase), ",
              "h_stride=$(m.horizontal_stride), ",
			  "h_phase=$(m.horizontal_phase), ",
			  "extensionAble=$(m.extensionAble))")
end

"""
    MoveSet

Defines the full set of allowed alignment moves used during pairwise or multiple sequence alignment.

# Fields
- `match_moves::Vector{Move}`: Moves that represent matches or substitutions between nucleotides. 
- `hor_moves::Vector{Move}`: Moves that introduce gaps in the **reference (horizontal)** sequence.
- `vert_moves::Vector{Move}`: Moves that introduce gaps in the **query (vertical)** sequence.

# Description
A `MoveSet` groups the allowable `Move`s into categories used during dynamic programming alignment. Each type controls how the algorithm can transition between states, including nucleotide-level match and gap moves with customizable reading frame behavior.

Used by alignment algorithms such as `seed_chain_align` and  `msa_codon_align` to control the scoring and allowed operations during alignment.
"""
struct MoveSet # TODO swap to tuples
	match_moves::Vector{Move}
	hor_moves::Vector{Move}
	vert_moves::Vector{Move}
end

MoveSet(; match_moves, hor_moves, vert_moves) = MoveSet(match_moves, hor_moves, vert_moves)

function show(io::IO, ms::MoveSet)
	print(io, "MoveSet(\n",
			  "  match_moves=$(ms.match_moves) \n",
			  "  hor_moves=$(ms.hor_moves) \n",
			  "  vert_moves=$(ms.vert_moves)\n)")
end

"""
    std_codon_moveset()

Return a standard codon-aware `MoveSet` for sequence alignment.

# Returns
- `MoveSet`: Contains codon-aware match and gap moves.

# Moves
- Match: `Move(1, 0.0)`
- Horizontal (gaps in reference sequence):
  - `Move(1, 2.0, 1, 0, 1, 0, false)`
  - `Move(3, 2.0, 1, 0, 3, 0, true)`
- Vertical (gaps in query sequence):
  - `Move(1, 2.0, 1, 0, 1, 0, false)`
  - `Move(3, 2.0, 1, 0, 3, 0, true)`

Use this as a default `MoveSet` for codon-preserving alignments.
"""
function std_codon_moveset()
	match_moves = [Move(1,.0)]
	hor_moves =  [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]
	vert_moves = [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,3,0, true)]
	return MoveSet(match_moves,hor_moves,vert_moves)
end
"""
    pairwise_noisy_moveset()

"""
function pairwise_noisy_moveset()
    match_moves = [Move(1,.0)]
    hor_moves =  [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,1,0, true)]
	vert_moves = [Move(1, 2.0, 1, 0, 1,0, false), Move(3, 2.0, 1,0,1,0, true)]
    return MoveSet(match_moves,hor_moves,vert_moves)
end

"""
    get_all_moves(moveset)

Collects the moves from the moveset. 
# Returns
    moveset.match_moves, moveset.hor_moves, moveset.vert_moves
"""

function get_all_moves(moveset::MoveSet)
	return moveset.match_moves, moveset.hor_moves, moveset.vert_moves
end