export Painter, paint,
    InitPainter,
    ObjectPainter,
    IDPainter

abstract type Painter end

function paint end

@with_kw struct InitPainter <: Painter
    path::String
    background::String = "#e7e7e7"
end

@with_kw struct ObjectPainter <: Painter
    dot_color = "#b4b4b4"
    # highlight = "#ea3433"
    probe_color = "#a0a0a0"
    wall_color = "black"
    alpha::Float64 = 1.0
end

function paint(p::ObjectPainter, dot::Dot)
    _draw_circle(get_pos(dot), dot.radius, p.dot_color,
                 opacity = p.alpha)
    return nothing
end

function paint(p::ObjectPainter, w::Wall)
    _draw_arrow(w.p1, w.p2, p.wall_color, arrowheadlength=0.0)
    return nothing
end

@with_kw struct IDPainter <: Painter
    label::Bool = true
    label_size::Float64 = 40.0
    alpha::Float64 = 1.0
end

function paint(p::IDPainter, d::Dot, v::Int64)
    pos = get_pos(d)
    p.label && _draw_text("$v", pos .+ [d.radius, d.radius],
                          size = p.label_size)
    return nothing
end
