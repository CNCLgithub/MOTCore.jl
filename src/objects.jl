export Thing,
    get_pos,
    Object,
    Dot,
    Wall,
    UniformEnsemble

abstract type Object end


"""
    get_pos(::Object)

The 2D position of an object.
"""
function get_pos end

"""
    get_vel(::Object)

The 2D instantaneous velocity of an object
"""
function get_vel end

"A point object - used for tracking targets"
@with_kw struct Dot <: Object
    # Dynamics
    radius::Float64
    # Kinematics
    pos::SVector{2, Float64}
    vel::SVector{2, Float64}
end

get_pos(d::Dot) = d.pos
get_vel(d::Dot) = d.vel
target(d::Dot) = d.target

"A wall (the boundry of the scene)"
struct Wall <: Object
    # Dynamics
    d::Float64 # the distance from the center
    normal::SVector{2, Float64} # normal vector
    nd::SVector{2, Float64}
    function Wall(d::Float64, normal::SVector{2, Float64})
        new(d, normal, d * normal)
    end
    # Kinematics <none>
    # Graphcis <none>
end

const WALL_ANGLES = [0.0, pi/2, pi, 3/2 * pi]

function init_walls(width::Float64)
   ws = Vector{Wall}(undef, 4)
   @inbounds for (i, theta) in enumerate(WALL_ANGLES)
       normal = SVector{2, Float64}([cos(theta), sin(theta)])
       sign = (-1)^(i > 2)
       ws[i] = Wall(width * 0.5 * sign, normal)
    end
    return SVector{4, Wall}(ws)
end
