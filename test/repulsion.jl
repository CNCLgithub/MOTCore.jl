using Gen
using JSON
using MOTCore
using Accessors: setproperties

function trial_to_dict(states::Vector)
    positions = map(states) do state
        map(state.objects) do dot
            get_pos(dot)
        end
    end
    Dict(:positions => positions)
end

function write_json(dataset::Vector, json_path::String)
    out = map(trial_to_dict, dataset)
    open(json_path, "w") do f
        write(f, json(out))
    end
    return nothing
end

function gen_trial(wm::RepulsionWM, vel_a::Float64, vel_b::Float64,
                   t::Int64, tw::Int64)
    split = uniform_discrete(-tw, tw)
    dur_a = t + split
    init, part_a =
        wm_repulsion(dur_a, setproperties(wm; vel = vel_a))
    dur_b = t - split
    part_b =
        repulsion_chain(dur_b,
                        last(part_a),
                        setproperties(wm; vel = vel_b))
    states = vcat(part_a, part_b)
end

function write_condlist(n::Int64, c::Int64, path::String)
    out = []
    for i = 1:c
        push!(out, [i-1, false])
        if i > (c - n)
            push!(out, [i-1, true])
        end
    end
    open(path, "w") do f
        write(f, json(out))
    end
    return nothing
end

function test()
    wm = RepulsionWM(;
                     n_dots = 8,
                     wall_repulsion = 200.0,
                     dot_repulsion = 200.0,
                     distance_factor = 15.0,
                     max_distance = 150.0,
                     )

    # t = 60
    t = 120
    tw = 48
    easy_vel = 6.5
    hard_vel = 9.0
    # n = 1
    n = 5
    c = 0
    dataset = []
    # easy only
    for i = 1:n
        trial = gen_trial(wm, easy_vel, easy_vel, t,tw)
        push!(dataset, trial)
        c += 1
    end

    # hard only
    for i = 1:n
        trial = gen_trial(wm, hard_vel, hard_vel, t,tw)
        push!(dataset, trial)
        c += 1
    end

    # easy-hard
    for i = 1:n
        trial = gen_trial(wm, easy_vel, hard_vel, t,tw)
        push!(dataset, trial)
        c += 1
    end
    write_json(dataset, "test/output/dataset.json")
    write_condlist(n, c, "test/output/trial_list.json")
    # write_json(dataset, "test/output/examples.json")
end

test();
