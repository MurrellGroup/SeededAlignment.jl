import Base: show

struct Move
    ref::Bool
    step_length::Int64
    score::Float64
    extendable::Bool

    function Move(ref::Bool,  step_length::Int64, score::Float64, extendable::Bool)
        (step_length in (1, 2, 3)) || throw(ArgumentError("step_length must be 1, 2, or 3"))
        #(score < 0) || throw(ArgumentError("score must be negative"))
        (!ref || (step_length == 3 && extendable)) || throw(ArgumentError("invalid ref:\n when ref=true it is required that step_length=3 and extendable=true."))
        new(ref, step_length, score, extendable)
    end
end
#TODO fix the useful constructor.
Move(; ref=false::Bool, step_length::Int64, score::Float64, extendable=false) = Move(ref, step_length, score, extendable)

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