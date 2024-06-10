export SchollWM,
    SchollState

@with_kw struct SchollWM <: WorldModel

    # EPISTEMICS
    n_dots::Int64 = 2
    n_targets::Int64 = ceil(Int64, n_dots * 0.5)
    target_p::Float64 = n_targets / n_dots

    # DYNAMICS
    dot_radius::Float64 = 20.0
    area_width::Float64 = 1200.0
    area_height::Float64 = 800.0
    dimensions::Tuple{Float64, Float64} = (area_width, area_height)
    vel::Float64 = 7.0 # base velocity
    "Amount to change velocity"
    vel_step::Float64 = 1.0
    "Probability of changing velocity"
    vel_prob::Float64 = 0.1
    vel_max::Float64 = 1.5 * vel
    vel_min::Float64 = 0.5 * vel
end

struct SchollState <: WorldState{SchollWM}
    objects::Vector{Dot}
    targets::Vector{Bool}
end

function SchollState(st::SchollState, dots::Array{<:Dot})
    setproperties(st; objects = dots)
end

function SchollState(gm::SchollWM, dots, targets)
    SchollState(dots, targets)
end

function tracker_bounds(gm::SchollWM)
    @unpack area_width, area_height, dot_radius = gm
    xs = (-0.5*area_width + dot_radius, 0.5*area_width - dot_radius)
    ys = (-0.5*area_height + dot_radius, 0.5*area_height - dot_radius)
    (xs, ys)
end

function Dot(wm::SchollWM, pos, vel)
    Dot(wm.dot_radius, pos, vel)
end

function step(gm::SchollWM,
              state::SchollState,
              deltas::AbstractVector{<:SVector{2, Float64}})
    # Dynamics (computing forces)
    @unpack (n_dots, area_height, area_width, dot_radius,
             vel_min, vel_max) = gm
    @unpack objects = state
    new_dots = Vector{Dot}(undef, n_dots)

    @inbounds for i in eachindex(objects)
        dot = objects[i]
        pos = get_pos(dot)
        # random delta for velocity
        new_vel = get_vel(dot) + deltas[i]
        # clamp velocity
        nv = normalize(new_vel)
        mag = clamp(norm(new_vel), vel_min, vel_max)
        new_vel = nv .* mag

        # adjust out of bounds positions
        new_pos = pos + new_vel
        if !in_bounds(gm, new_pos)
            new_vel = deflect(gm, pos, new_vel)
            new_pos = pos + new_vel
        end
        new_pos = SVector{2, Float64}(
            clamp(new_pos[1], -area_width * 0.5 + dot_radius,
                  area_width * 0.5  - dot_radius),
            clamp(new_pos[2], -area_height * 0.5 + dot_radius,
                  area_height * 0.5  - dot_radius)
        )

        new_dots[i] = setproperties(dot; pos = new_pos, vel = new_vel)
    end
    SchollState(state, new_dots)
end

function step(gm::SchollWM,
              state::SchollState)
    n_dots = gm.n_dots
    step(gm, state, Fill(zero2vec, n_dots))
end


include("gen.jl")

################################################################################
# Misc
################################################################################

function in_bounds(gm::SchollWM, position::SVector{2, Float64})
    @unpack area_height, area_width, dot_radius = gm
    (abs(position[1]) + dot_radius) < 0.5 * area_width &&
        (abs(position[2]) + dot_radius) < 0.5 * area_height
end

function deflect(gm::SchollWM, pos::T, vel::T) where {T<:SVector{2, Float64}}
    @unpack area_height, area_width, dot_radius = gm
    tpos = pos + vel
    nv = normalize(vel)
    dir = if tpos[1] < (-0.5 * area_width + dot_radius)
        # left wall
        vec_down
    elseif tpos[1] > (0.5 * area_width - dot_radius)
        # right wall
        vec_up
    elseif tpos[2] < (-0.5 * area_height + dot_radius)
        # bottom wall
        vec_right
    else
        # top wall
        vec_left
    end
    aoi = acos(dot(nv, dir))
    new_vel = rot2dvec(vel, 2.0 * aoi)
end

function objects_from_positions(gm::SchollWM, positions)
    nx = length(positions)
    dots = Vector{Dot}(undef, nx)
    for i = 1:nx
        xy = Float64.(positions[i][1:2])
        dots[i] = Dot(gm,
                      SVector{2}(xy),
                      SVector{2, Float64}(zeros(2)))
    end
    return dots
end

function state_from_positions(gm::SchollWM, positions, targets)
    nt = length(positions)
    states = Vector{SchollState}(undef, nt)
    for t = 1:nt
        if t == 1
            dots = objects_from_positions(gm, positions[t])
            states[t] = SchollState(gm, dots, targets)
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
            new_dots[i] = setproperties(objects[i]; pos = new_pos, vel = new_vel)
        end
        states[t] = SchollState(gm, new_dots, targets)
    end
    return states
end

function paint(p::InitPainter, wm::SchollWM, st::SchollState)
    @unpack area_width, area_height = wm
    Drawing(area_width, area_height, p.path)
    Luxor.origin()
    background(p.background)
end

function paint(p::Painter, st::SchollState)
    for o in st.objects
        paint(p, o)
    end
    return nothing
end

# function paint(p::Union{IDPainter,KinPainter},
function paint(p::IDPainter,
               st::SchollState)
    for i in eachindex(st.objects)
        paint(p, st.objects[i], i)
    end
    return nothing
end
