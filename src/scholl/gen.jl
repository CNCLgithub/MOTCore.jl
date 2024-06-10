export wm_scholl,
    scholl_chain,
    scholl_init,
    scholl_kernel

################################################################################
# Initial State
################################################################################
@gen static function scholl_dot(wm::SchollWM)
    xs, ys = tracker_bounds(wm)
    x = @trace(uniform(xs[1], xs[2]), :x)
    y = @trace(uniform(ys[1], ys[2]), :y)

    ang = @trace(uniform(0.0, 2*pi), :ang)
    mag = @trace(normal(wm.vel, wm.vel_step), :std)

    pos = SVector{2, Float64}(x, y)
    vel = SVector{2, Float64}(mag*cos(ang), mag*sin(ang))

    new_dot::Dot = Dot(wm, pos, vel)
    return new_dot
end

@gen static function target_flip(w::Float64)
    f::Bool = @trace(bernoulli(w), :target)
    return f
end

@gen static function scholl_init(wm::SchollWM)
    wms = Fill(wm, wm.n_dots)
    dots = @trace(Gen.Map(scholl_dot)(wms), :dots)
    target_ps = Fill(wm.target_p, wm.n_dots)
    targets = @trace(Gen.Map(target_flip)(target_ps), :targets)
    state::SchollState = SchollState(wm, dots, targets)
    return state
end

################################################################################
# Dynamics
################################################################################


@gen (static) function scholl_delta(wm::SchollWM)
    flip = @trace(bernoulli(wm.vel_prob), :flip)
    fx = @trace(normal(0, wm.vel_step), :fx)
    fy = @trace(normal(0, wm.vel_step), :fy)
    f::SVector{2, Float64} = SVector{2, Float64}(fx * flip,
                                                 fy * flip)
    return f
end

@gen (static) function scholl_kernel(t::Int,
                                        prev_st::SchollState,
                                        wm::SchollWM)
    deltas = @trace(Gen.Map(scholl_delta)(Fill(wm, wm.n_dots)), :deltas)
    next_st::SchollState = step(wm, prev_st, deltas)
    return next_st
end

const scholl_chain = Gen.Unfold(scholl_kernel)


@gen (static) function wm_scholl(k::Int, wm::SchollWM)
    init_state = @trace(scholl_init(wm), :init_state)
    states = @trace(scholl_chain(k, init_state, wm), :kernel)
    result = (init_state, states)
    return result
end
