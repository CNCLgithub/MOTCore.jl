using Gen
using JSON
using MOTCore
using Statistics: mean
using LinearAlgebra: norm
using Accessors: setproperties

function trial_to_dict(states::Vector)
    positions = map(states) do state
        map(state.objects) do dot
            get_pos(dot)
        end
    end
    Dict(:positions => positions)
end

function write_dataset(dataset::Vector, json_path::String)
    out = map(trial_to_dict, dataset)
    open(json_path, "w") do f
        write(f, json(out))
    end
    return nothing
end

function write_condlist(condlist::Vector, path::String)
    open(path, "w") do f
        write(f, json(condlist))
    end
    return nothing
end

function gen_init(wm::RepulsionWM, dur::Int64, low_tdd, high_tdd)
    accept = false
    local part
    local tdd::Float64
    while !accept
        state = burnin(wm, 12)
        part = repulsion_chain(dur, state, wm)
        tdd = mean(tddensity, part)
        accept = tdd > low_tdd && tdd < high_tdd
        print("tdd: $(tdd) \r")
    end
    return part, tdd
end

function gen_part(st::RepulsionState, wm::RepulsionWM, dur::Int64, low_tdd, high_tdd)
    accept = false
    local part
    local tdd::Float64
    while !accept
        part = repulsion_chain(dur, st, wm)
        tdd = mean(tddensity, part)
        accept = tdd > low_tdd && tdd < high_tdd
        print("tdd: $(tdd) \r")
    end
    return part, tdd
end

function gen_trial(wm::RepulsionWM, tot_dur::Int64, parts)
    nparts = length(parts)
    dur_part = Int64(tot_dur / nparts)
    # burnin period to prevent objects from occluding
    trial = RepulsionState[]
    tdds = Vector{Float64}(undef, nparts)
    for (i, (low, high)) = enumerate(parts)
        local part
        local tdd
        @time if i == 1
            (part, tdd) = gen_init(wm, dur_part, low, high)
        else
            (part, tdd) = gen_part(last(trial), wm, dur_part, low, high)
        end
        trial = vcat(trial, part)
        tdds[i] = tdd
    end
    return trial, tdds
end

function tddensity(state::RepulsionState)
    # first 4 are targets
    objects = state.objects
    # avg_tdd = 0.0
    avg_tdd = Inf
    @inbounds for i = 1:4
        # tdd = Inf
        tpos = get_pos(objects[i])
        for j = 5:8
            d = norm(tpos - get_pos(objects[j]))
            # avg_tdd += d # min(tdd, d)
            avg_tdd = min(avg_tdd, d)
        end
        # avg_tdd += tdd
    end
    # avg_tdd / 16
    avg_tdd
end

function burnin(wm::RepulsionWM, steps::Int64)
    state = repulsion_init(wm)
    for t = 1:steps
        state = repulsion_kernel(t, state, wm)
    end
    return state
end

function warmup(wm::RepulsionWM, k::Int64=240, n::Int64=1000)
    m = tddensity(burnin(wm, 12))
    v = 0.0
    c = 1
    for _ = 1:n
        state = burnin(wm, 12)
        part = repulsion_chain(k, state, wm)
        x = mean(tddensity, part)
        m_prev = m
        m = m_prev + ((x - m_prev) / c)
        v += (x - m_prev) * (x - m)
        c += 1
    end
    (m, sqrt(v / (c - 1)))
end

function test()
    wm = RepulsionWM(;
        n_dots=8,
        wall_repulsion=10.0,
        dot_repulsion=12.0,
        distance_factor=37.0,
        max_distance=200.0,
        vel=7.0,
        max_accel=4.0,
        max_vel=10.0,
        force_sd=1.0,
    )

    # dataset parameters
    tot_dur = 360  # 15s

    @time tdd = tdd_mu, tdd_sd = warmup(wm, tot_dur, 1000)
    D = [0.0, tdd_mu - 4.00 * tdd_sd]
    M = [tdd_mu - 1.00 * tdd_sd, tdd_mu + 1.00 * tdd_sd]
    E = [tdd_mu + 6.00 * tdd_sd, Inf]

    @show tdd
    @show D
    @show M
    @show E

    dataset = []
    cond_list = []
    base = "test/output/dataset"
    isdir(base) || mkdir(base)

    # trial = gen_trial(wm, tot_dur, [[0.0, Inf]])
    # render_scene(wm, trial, "$(base)/test")

    perms = [[D, D, D],
             [E, E, E],
             [D, M, E],
             [D, M, E],
             [E, M, D],
             [E, M, D]]

    df = DataFrame(:scene => Int32[],
                   :epoch => Int32[],
                   :tdd => Float32[])
    # total of 12 trials
    for (i, perm) = enumerate(perms)
        ne = length(perm)
        trial, epochs = gen_trial(wm, tot_dur, perm)
        @show perm
        @show epochs
        push!(dataset, trial)
        push!(cond_list, [i, false])
        push!(cond_list, [i, true])
        append!(df, Dict(:scene => i,
                         :epoch => 1:ne,
                         :tdd => epochs))
    end
    write_dataset(dataset, "test/output/dataset.json")
    write_condlist(cond_list, "test/output/trial_list.json")
    CSV.write("test/output/tdd.csv", df)
end

test();
