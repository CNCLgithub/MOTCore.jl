export wm_repulsion,
    repulsion_chain

################################################################################
# Initial State
################################################################################
@gen static function repulsion_dot(wm::RepulsionWM)
    xs, ys = tracker_bounds(wm)
    x = @trace(uniform(xs[1], xs[2]), :x)
    y = @trace(uniform(ys[1], ys[2]), :y)

    ang = @trace(uniform(0.0, 2*pi), :ang)
    mag = @trace(normal(wm.vel, 1e-2), :std)

    pos = SVector{2, Float64}([x, y])
    vel = SVector{2, Float64}([mag*cos(ang), mag*sin(ang)])

    new_dot::Dot = Dot(wm, pos, vel)
    return new_dot
end

@gen static function target_flip(w::Float64)
    f::Bool = @trace(bernoulli(w), :target)
    return f
end

@gen static function repulsion_init(wm::RepulsionWM)
    wms = Fill(wm, wm.n_dots)
    dots = @trace(Gen.Map(repulsion_dot)(wms), :dots)
    target_ps = Fill(wm.target_p, wm.n_dots)
    targets = @trace(Gen.Map(target_flip)(target_ps), :targets)
    state::RepulsionState = RepulsionState(wm, dots, targets)
    return state
end

################################################################################
# Dynamics
################################################################################

@gen (static) function repulsion_kernel(t::Int,
                                        prev_st::RepulsionState,
                                        wm::RepulsionWM)
    next_st::RepulsionState = step(wm, prev_st)
    return next_st
end

const repulsion_chain = Gen.Unfold(repulsion_kernel)


@gen (static) function wm_repulsion(k::Int, wm::RepulsionWM)
    init_state = @trace(repulsion_init(wm), :init_state)
    states = @trace(repulsion_chain(k, init_state, wm), :kernel)
    result = (init_state, states)
    return result
end
