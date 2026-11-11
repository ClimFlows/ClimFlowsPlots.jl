module SpectralSphere

"""
    orthographic(lons, lats, field; options...)

Plot a longitude–latitude field on an orthographic globe. Coordinates are in
degrees and longitude indexes the first dimension of `field`. Pass a Makie
`Observable` to update the field in place during animation. Load Makie and
GeoMakie to enable this method. Keyword arguments are passed to `surface!`.
"""
function orthographic end

# convenience functions

@views function upscale(field)
    field = 0.5 * (field[1:2:end, :] + field[2:2:end, :])
    field = 0.5 * (field[:, 1:2:end] + field[:, 2:2:end])
end

bounds_lon(lons) = bounds(lons[end] - 360, lons, lons[1] + 360)
bounds_lat(lats) = bounds(-90.0, lats, 90.0)
function bounds(x0, x, xend)
    x = (x[2:end] + x[1:end-1]) / 2
    return [x0, x..., xend]
end

end
