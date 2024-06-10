
const zero2vec = SVector{2, Float64}(0., 0.)

function rot2dvec(vec::SVector{2, Float64}, rad::Float64)
    cs = cos(rad)
    sn = sin(rad)
    x, y = vec
    SVector{2, Float64}(x * cs - y * sn,
                        x * sn + y * cs)
end

const vec_up = SVector{2, Float64}([0, 1])
const vec_down = SVector{2, Float64}([0, -1])
const vec_left = SVector{2, Float64}([-1, 0])
const vec_right = SVector{2, Float64}([1, 0])
