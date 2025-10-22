module VoronoiSphere

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

end # VoronoiSphere
