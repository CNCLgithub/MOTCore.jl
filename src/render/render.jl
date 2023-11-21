using Luxor

export render_scene

include("helpers.jl")
include("painters/painters.jl")

function render_scene(wm::T,
                      states::AbstractArray{<:WorldState{<:T}},
                      base::String) where {T<:WorldModel}

    isdir(base) && rm(base, recursive=true)
    mkdir(base)

    nt = length(states)

    objp = ObjectPainter()
    idp = IDPainter()

    @inbounds for i = 1:nt
        print("rendering scene... timestep $i / $nt \r")
        init = InitPainter(path = "$base/$i.png",
                           background = "white")
        state = states[i]
        paint(init, wm, state)
        paint(objp, state)
        paint(idp, state)
        finish()
    end
    return nothing
end
