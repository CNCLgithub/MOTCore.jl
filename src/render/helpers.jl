"""
    helper to draw text
"""
function _draw_text(text, position; opacity=1.0, color="black", size=30)
    setopacity(opacity)
    sethue(color)
    Luxor.fontsize(size)
    point = Luxor.Point(position[1], -position[2])
    Luxor.text(text, point, halign=:right, valign=:bottom)
end

"""
    helper to draw circle
"""
function _draw_circle(position, radius, color;
                      opacity=1.0, style=:fill,
                      pattern="solid", line=5)
    if style==:stroke
        setline(line)
        setdash(pattern)
    end
    setopacity(opacity)
    sethue(color)
    point = Luxor.Point(position[1], -position[2])
    Luxor.circle(point, radius, style)
end

"""
    helper to draw array (used to draw the predicted tracker masks)
"""
function _draw_array(array, area_width, area_height, grid_width, grid_height, color; opacity=1.0)
    sethue(color)

    #tiles = Tiler(gm.area_width, gm.area_height, gm.img_width, gm.img_height, margin=0)
    tiles = Tiler(area_height, area_width, grid_height, grid_width, margin=0)

    for (pos, n) in tiles
        # reading value from the array
        row = tiles.currentrow
        col = tiles.currentcol
        value = array[row, col]

        # scaling opacity according to the value
        setopacity(opacity*value)

        box(pos, tiles.tilewidth, tiles.tileheight, :fill)
    end
end

"""
    helper to draw arrow
"""
function _draw_arrow(startpoint, endpoint, color;
                     opacity=1.0,
                     linewidth=5.0, arrowheadlength=15.0)
    setopacity(opacity)
    sethue(color)
    p1 = Luxor.Point(startpoint[1], -startpoint[2])
    p2 = Luxor.Point(endpoint[1], -endpoint[2])
    Luxor.arrow(p1, p2, linewidth=linewidth, arrowheadlength=arrowheadlength)
end

function _draw_line(startpoint, endpoint, color;
                    opacity=1.0,
                    linewidth=5.0)
    setopacity(opacity)
    sethue(color)
    p1 = Luxor.Point(startpoint[1], -startpoint[2])
    p2 = Luxor.Point(endpoint[1], -endpoint[2])
    Luxor.line(p1, p2, :stroke)
end
