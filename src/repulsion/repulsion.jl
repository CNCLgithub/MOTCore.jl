export RepulsionWM,
    RepulsionState

@with_kw struct RepulsionWM <: WorldModel

    # EPISTEMICS
    n_dots::Int64 = 2
    n_targets::Int64 = ceil(Int64, n_dots * 0.5)
    target_p::Float64 = n_targets / n_dots

    # DYNAMICS
    dot_radius::Float64 = 20.0
    area_width::Float64 = 800.0
    area_height::Float64 = 800.0
    dimensions::Tuple{Float64, Float64} = (area_width, area_height)
    dot_repulsion::Float64 = 80.0
    wall_repulsion::Float64 = 50.0
    distance_factor::Float64 = 10.0
    max_distance::Float64 = 100.0
    vel::Float64 = 10.0 # base velocity
    rep_inertia::Float64 = 0.9
end

struct RepulsionState <: WorldState{RepulsionWM}
    walls::SVector{4, Wall}
    objects::Vector{Dot}
    targets::Vector{Bool}
end

function RepulsionState(st::RepulsionState, dots::Array{<:Dot})
    setproperties(st; objects = dots)
end

function RepulsionState(gm::RepulsionWM, dots, targets)
    walls = init_walls(gm.area_width)
    RepulsionState(walls, dots, targets)
end

function tracker_bounds(gm::RepulsionWM)
    @unpack area_width, area_height, dot_radius = gm
    xs = (-0.5*area_width + dot_radius, 0.5*area_width - dot_radius)
    ys = (-0.5*area_height + dot_radius, 0.5*area_height - dot_radius)
    (xs, ys)
end

function Dot(wm::RepulsionWM, pos, vel)
    Dot(wm.dot_radius, pos, vel)
end

function step(gm::RepulsionWM,
              state::RepulsionState)

    # Dynamics (computing forces)
    n_dots = gm.n_dots
    @unpack walls, objects = state
    new_dots = Vector{Dot}(undef, n_dots)

    @inbounds for i in eachindex(objects)
        dot = objects[i]
        # force accumalator
        facc = MVector{2, Float64}(zeros(2))
        # interactions with walls
        for w in walls
            force!(facc, gm, w, dot)
        end
        for j in eachindex(objects)
            i === j && continue
            force!(facc, gm, objects[j], dot)
        end
        new_dots[i] = update_kinematics(gm, dot, facc)
    end
    RepulsionState(state, new_dots)
end

function force!(f::MVector{2, Float64}, gm::RepulsionWM, w::Wall, d::Dot)
    pos = get_pos(d)
    @unpack wall_repulsion, max_distance, distance_factor = gm
    v = LinearAlgebra.norm(w.normal .* pos + w.nd)
    v > max_distance && return nothing
    mag = wall_repulsion * exp(-v/distance_factor)
    f .+= mag * w.normal
    return nothing
end

function force!(f::MVector{2, Float64}, gm::RepulsionWM, x::Dot, d::Dot)
    @unpack dot_repulsion, max_distance, distance_factor = gm
    v = get_pos(d) - get_pos(x)
    nv = norm(v)
    nv > max_distance && return nothing
    mag = dot_repulsion * exp(-nv/distance_factor)
    delta_f = mag * (v./nv)
    f .+= delta_f
    return nothing
end

function update_kinematics(gm::RepulsionWM, d::Dot, f::MVector{2, Float64})
    @unpack rep_inertia, vel, area_height = gm
    new_vel = (d.vel ./ norm(d.vel)) * vel
    new_vel += (rep_inertia * f)
    new_pos = clamp.(get_pos(d) + new_vel,
                     -area_height * 0.5 + d.radius,
                     area_height * 0.5  - d.radius)
    any(isnan, new_pos) && error()
    setproperties(d; pos = new_pos, vel = new_vel)
end

include("gen.jl")

################################################################################
# Misc
################################################################################

function objects_from_positions(gm::RepulsionWM, positions, targets)
    nx = length(positions)
    dots = Vector{Dot}(undef, nx)
    for i = 1:nx
        xy = Float64.(positions[i][1:2])
        dots[i] = Dot(gm,
                      SVector{2}(xy),
                      SVector{2, Float64}(zeros(2)),
                      Bool(targets[i]))
    end
    return dots
end

function state_from_positions(gm::RepulsionWM, positions, targets)
    nt = length(positions)
    states = Vector{RepulsionState}(undef, nt)
    for t = 1:nt
        if t == 1
            dots = objects_from_positions(gm, positions[t], targets)
            states[t] = RepulsionState(gm, dots)
            continue
        end

        prev_state = states[t-1]
        @unpack objects = prev_state
        ni = length(objects)
        new_dots = Vector{Dot}(undef, ni)
        for i = 1:ni
            old_pos = Float64.(positions[t-1][i][1:2])
            new_pos = SVector{2}(Float64.(positions[t][i][1:2]))
            new_vel = SVector{2}(new_pos .- old_pos)
            new_dots[i] = sync_update(objects[i], KinematicsUpdate(new_pos, new_vel))
        end
        states[t] = RepulsionState(gm, new_dots)
    end
    return states
end

function paint(p::InitPainter, wm::RepulsionWM, st::RepulsionState)
    height, width = wm.dimensions
    Drawing(width, height, p.path)
    Luxor.origin()
    background(p.background)
end

function paint(p::Painter, st::RepulsionState)
    for o in st.objects
        paint(p, o)
    end
    return nothing
end

# function paint(p::Union{IDPainter,KinPainter},
function paint(p::IDPainter,
               st::RepulsionState)
    for i in eachindex(st.objects)
        paint(p, st.objects[i], i)
    end
    return nothing
end
