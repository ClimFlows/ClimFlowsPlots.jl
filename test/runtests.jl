using ClimFlowsPlots
using ClimFlowsPlots: VoronoiSphere as VSPlots
using Test

using ClimFlowsData: DYNAMICO_reader, DYNAMICO_meshfile
using VoronoiSpheres: VoronoiSphere
using GeoMakie, CairoMakie
using NetCDF: ncread

fun(lon, lat) = cos(lat)^4*cos(4*lon) # in radians
savefig(name, fig) = save(joinpath(savepath, name), fig)

savepath = mkpath(joinpath(@__DIR__, "artifacts"))
meshname = "uni.2deg.mesh.nc"
sphere = VoronoiSphere(DYNAMICO_reader(ncread, DYNAMICO_meshfile(meshname)) ; prec=Float32)
pv = GeoMakie.Observable(fun.(sphere.lon_i, sphere.lat_i))

@testset "ClimFlowsPlots.jl" begin
    savefig("orthographic.png", VSPlots.plot_orthographic(sphere, pv ; colormap=:berlin))
    @test true
    savefig("2D.png", VSPlots.plot_2D(sphere, pv ; colormap=:berlin))
    @test true
    savefig("native_3D.png", VSPlots.plot_native_3D(sphere, pv ; colormap=:berlin))
    @test true
end
