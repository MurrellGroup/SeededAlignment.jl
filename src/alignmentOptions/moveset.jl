import Base: show

struct Move
    ref::Bool
    step_length::Int64
    score::Float64
    extendable::Bool

    function Move(ref::Bool,  step_length::Int64, score::Float64, extendable::Bool)
        (step_length in (1, 2, 3)) || throw(ArgumentError("step_length must be 1, 2, or 3"))
        #(score < 0) || throw(ArgumentError("score must be negative"))
        (!ref || (step_length == 3 && extendable)) || throw(ArgumentError("Invalid ref move:\n when ref=true it is required that step_length=3 and extendable=true."))
        new(ref, step_length, score, extendable)
    end
end
# useful constructors
Move(; ref=false::Bool, step_length::Int64, score::Float64, extendable=false) = Move(ref, step_length, score, extendable)
# assumes ref true
Move(; score::Float64) = Move(ref=true, step_length=3, score=score, extendable=true)

# stack allocated 
struct Moveset{X,Y}
    # gap operation in top/first provided sequence
    vert_moves::NTuple{X,Move}
    # gap operations in left/second provided sequence
    hor_moves::NTuple{Y,Move}
end
# reference informed Moveset constructor
Moveset(; ref_insertions, ref_deletions) = Moveset(ref_insertions, ref_deletions)
# non-reference informed Moveset constructor - alignment operations are symmetric
Moveset(gap_moves) = Moveset(gap_moves, gap_moves)

function show(io::IO, ms::Moveset)
	println(io, "Moveset(\n","vert_moves=$(ms.vert_moves)\n)","hor_moves=$(ms.hor_moves)\n)")
end

function contains_ref_move(moveset::Moveset)
    ref_move_found = any(move.ref for move in moveset.vert_moves) || any(move.ref for move in moveset.hor_moves)
    return ref_move_found
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
    vert_moves = (Move(ref=false, step_length=1, score=2.0, extendable=false), Move(ref=true, step_length = 3, score=2.0, extendable=true))
    hor_moves =  (Move(ref=false, step_length=1, score=2.0, extendable=false), Move(ref=true, step_length = 3, score=2.0, extendable=true))
    return Moveset(vert_moves,hor_moves)
end
"""
    pairwise_noisy_moveset()

"""
function std_noisy_moveset()
    vert_moves = (Move(ref=false, step_length=1, score=2.0, extendable=false), Move(ref=false, step_length = 3, score=2.0, extendable=true))
    hor_moves =  (Move(ref=false, step_length=1, score=2.0, extendable=false), Move(ref=false, step_length = 3, score=2.0, extendable=true))
    return Moveset(vert_moves,hor_moves)
end