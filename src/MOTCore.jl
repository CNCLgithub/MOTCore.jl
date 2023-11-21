module MOTCore

export WorldModel,
    WorldState,
    step


using Gen
using Accessors
using Parameters
using FillArrays
using StaticArrays
using LinearAlgebra


abstract type WorldModel end

abstract type WorldState{W<:WorldModel} end

function step end

include("utils.jl")
include("objects.jl")
include("render/render.jl")
include("repulsion/repulsion.jl")

end # module MOTCore
