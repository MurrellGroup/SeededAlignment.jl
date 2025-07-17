import Base: show

struct Move
    ref::Bool
    step_length::Int64
    score::Float64
    extendable::Bool
    # force construction invariants
    function Move(ref::Bool,  step_length::Int64, score::Float64, extendable::Bool)
        (step_length in (1, 2, 3)) || throw(ArgumentError("step_length must be 1, 2, or 3"))
        #(score < 0) || throw(ArgumentError("score must be negative"))
        (!ref || (step_length == 3 && extendable)) || throw(ArgumentError("Invalid ref move:\n when ref=true it is required that step_length=3 and extendable=true."))
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

# stack allocated 
struct Moveset{X,Y}
    # gap operation in top/first provided sequence
    vert_moves::NTuple{X,Move}
    # gap operations in left/second provided sequence
    hor_moves::NTuple{Y,Move}
end
# reference informed Moveset constructor
function Moveset(; ref_insertions, ref_deletions) 
    Moveset(ref_insertions, ref_deletions)
end
# non-reference informed Moveset constructor - alignment operations are symmetric
Moveset(gap_moves) = Moveset(gap_moves, gap_moves)

function show(io::IO, ms::Moveset)
	println(io, "Moveset(\n","vert_moves=$(ms.vert_moves)\n)","hor_moves=$(ms.hor_moves)\n)")
end

function contains_ref_move(moveset::Moveset)
    ref_move_found = any(move.ref for move in moveset.vert_moves) || any(move.ref for move in moveset.hor_moves)
    return ref_move_found
end

# default movesets
# deletions same as insertions
const STD_CODON_MOVESET = Moveset(
    (
        Move(ref=false, step_length=1, score=2.0, extendable=false),
        Move(ref=true,  step_length=3, score=2.0, extendable=true)
    )
)
# deletions same as insertions
const STD_NOISY_MOVESET = Moveset(
    (
        Move(ref=false, step_length=1, score=2.0, extendable=false),
        Move(ref=false, step_length=3, score=2.0, extendable=true)
    )
)