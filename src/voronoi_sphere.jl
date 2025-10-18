module VoronoiSphere

# import GeometryBasics as GB
# import Makie
# using GeoMakie
# using ColorSchemes

# import ..SphericalInterpolations as SI

const DATA_ON_PRIMAL = "Data must be a scalar defined at primal cells of a spherical Voronoi mesh and passed as a Makie.Observable"

"""
    data = Makie.Observable(data) # must be on primal mesh
    fig = plot_2D(sphere, data ; resolution=1.0)
    display(fig)
    data[] = new_data # updates the plot
    display(fig)

Plots scalar data as a function of longitude and latitude. $DATA_ON_PRIMAL. Data is interpolated linearly to
a regular lon-lat grid of resolution `resolution` before plotting.
"""
function plot2D end # implemented in extension

"""
    data = Makie.Observable(data) # must be on primal mesh
    fig = plot_voronoi_orthographic(sphere, data ; zoom=1.6)
    display(fig)
    data[] = new_data # updates the plot
    display(fig)

Plots scalar data in a stereographic view. $DATA_ON_PRIMAL. Data is interpolated linearly to
a regular lon-lat grid of resolution `resolution` before plotting.
"""
function plot_orthographic end # implemented in extension

"""
    data = Makie.Observable(data) # must be on primal mesh
    fig = plot_native_3D(sphere, data ; zoom=1.6)
    display(fig)
    data[] = new_data # updates the plot
    display(fig)

    Plots scalar data in a stereographic view. $DATA_ON_PRIMAL. Data on
    native mesh is used without interpolation, but plotting interpolates linearly between
    primal cell centers (which coincides with triangle vertices).
"""
function plot_native_3D end # implemented in extension

#=
function plot_2D(sphere, data::Makie.Observable; resolution = 1.0, options...)
    # check that data is on primal mesh and interpolate to lon-lat
    @assert length(data[]) == length(sphere.Ai)
    lons, lats = -180:resolution:180, -90:resolution:90
    color = map(transpose ∘ SI.lonlat_interp(sphere, lons, lats), data)
    # create and return plot
    fig = Makie.Figure(size = (1440, 720))
    ax = Makie.Axis(fig[1, 1], aspect = 2)
    Makie.heatmap!(ax, color; colormap = :hot, options...)
    return fig
end

function plot_orthographic(sphere, data::Makie.Observable; resolution = 0.5, options...)
    # check that data is on primal mesh and interpolate to lon-lat
    @assert length(data[]) == length(sphere.Ai)
    lons, lats = -180:resolution:180, -90:resolution:90
    color = map(transpose ∘ SI.lonlat_interp(sphere, lons, lats), data)
    # create and return plot
    fig = Figure()
    ga = GeoAxis(fig[1, 1]; dest = "+proj=ortho +lon_0=19 +lat_0=50")
    surface!(ga, lons, lats, color; shading = NoShading, options...)
    lines!(ga, GeoMakie.coastlines())
    return fig
end

function plot_native_3D(sphere, data::Makie.Observable; zoom = 1.6, options...)
    # check that data is on primal mesh
    @assert length(data[]) == length(sphere.Ai)
    # build graphical mesh
    lon, lat, vertex = sphere.lon_i, sphere.lat_i, sphere.dual_vertex
    xyz(lon, lat) = cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)
    nodes = [GB.Point3f(xyz(lon[i], lat[i])) for i in eachindex(lon)]
    faces = [
        GB.GLTriangleFace((vertex[1, i], vertex[2, i], vertex[3, i])) for
        i in axes(vertex, 2)
    ]
    makiemesh = GB.Mesh(nodes, faces)
    # create and return plot
    fig, ax, obj = Makie.mesh(makiemesh; color = data, options...)
    Makie.scale!(ax.scene, zoom, zoom, zoom)
    return fig
end
=#

end # VoronoiSphere
