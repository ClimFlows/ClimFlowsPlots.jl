module SpectralSphere

# using GeoMakie

# orthographic plot

function orthographic(lons, lats, field; options...)
    fig = Figure()
    ga = GeoAxis(fig[1, 1]; dest = "+proj=ortho +lon_0=19 +lat_0=50")
    surface!(ga, lons, lats, field; shading = NoShading, options...)
    lines!(ga, GeoMakie.coastlines())
    return fig
end

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
