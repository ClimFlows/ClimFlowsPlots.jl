module GeoMakie_Ext

import GeometryBasics as GB
using Makie: Makie, lines!, surface!, NoShading
import GeoMakie
import ColorSchemes

import ClimFlowsPlots.VoronoiSphere: plot_2D, plot_orthographic, plot_native_3D
import ClimFlowsPlots.SphericalInterpolations as SI

function plot_2D(sphere, data::Makie.Observable, tree=SI.spherical_tree(sphere); resolution = 1.0, options...)
    # check that data is on primal mesh and interpolate to lon-lat
    @assert length(data[]) == length(sphere.Ai)
    lons, lats = -180:resolution:180, -90:resolution:90
    color = map(transpose ∘ SI.lonlat_interp(lons, lats, sphere, tree), data)
    # create and return plot
    fig = Makie.Figure(size = (1440, 720))
    ax = Makie.Axis(fig[1, 1], aspect = 2)
    Makie.heatmap!(ax, color; colormap = :hot, options...)
    return fig
end

function plot_orthographic(sphere, data::Makie.Observable, tree=SI.spherical_tree(sphere); resolution = 0.5, options...)
    # check that data is on primal mesh and interpolate to lon-lat
    @assert length(data[]) == length(sphere.Ai)
    lons, lats = -180:resolution:180, -90:resolution:90
    color = map(transpose ∘ SI.lonlat_interp(lons, lats, sphere, tree), data)
    # create and return plot
    fig = Makie.Figure()
    ga = GeoMakie.GeoAxis(fig[1, 1]; dest = "+proj=ortho +lon_0=19 +lat_0=50")
    surface!(ga, lons, lats, color; shading = NoShading, options...)
    lines!(ga, GeoMakie.coastlines())
    return fig
end

function plot_native_3D(sphere, data::Makie.Observable; zoom = 1.6, options...)
    # check that data is on primal mesh
    @assert length(data[]) == length(sphere.Ai)
    # build graphical mesh
    nodes = map(xyz, sphere.lon_i, sphere.lat_i)
    makiemesh = GB.Mesh(nodes, faces(sphere.dual_vertex))
    # create and return plot
    fig, ax, _ = Makie.mesh(makiemesh; color = data, options...)
    Makie.scale!(ax.scene, zoom, zoom, zoom)
    return fig
end
xyz(lon, lat) = GB.Point3f(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
faces(vertex::Matrix) = [ GB.GLTriangleFace((vertex[1, i], vertex[2, i], vertex[3, i])) for i in axes(vertex, 2) ] 
faces(vertex::Vector) = map(GB.GLTriangleFace, vertex)

end # module
