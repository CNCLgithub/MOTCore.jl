using MOTCore
using MOTCore: state_from_positions, render_scene
using JSON

wm = SchollWM(;
              n_dots=8,
              area_width = 720.0,
              area_height = 480.0,
              dot_radius = 20.0,
              vel=3.0,
              vel_step = 0.75,
              vel_prob = 0.15
)

targets = [true, true, true, true, false, false, false, false]


function main()
    dname = "0.5.1"
    local dataset
    open("test/output/$(dname)/dataset.json", "r") do f
        dataset = JSON.parse(f)
    end

    outpath = "test/output/$(dname)/dataset_renders"
    isdir(outpath) || mkdir(outpath)

    ntrials = length(dataset)
    for i = 1:ntrials
        positions = dataset[i]["positions"]
        states = state_from_positions(wm, positions, targets)
        scene_path = "$(outpath)/$(i)"
        isdir(scene_path) || mkdir(scene_path)
        render_scene(wm, states, scene_path)
    end

    return nothing
end

main();
