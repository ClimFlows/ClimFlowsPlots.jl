using ClimFlowsPlots
using ClimFlowsPlots: VoronoiSphere as VSPlots, SpectralSphere as SSPlots
using Test

using ClimFlowsData: DYNAMICO_reader, VoronoiMeshes
using VoronoiSpheres: VoronoiSphere
using GeoMakie, CairoMakie
using NetCDF: ncread


fun(lon, lat) = cos(lat)^4*cos(4*lon) # in radians
savefig(name, fig) = save(joinpath(savepath, name), fig)

savepath = mkpath(joinpath(@__DIR__, "artifacts"))
meshname = "uni.2deg.mesh.nc"
sphere = VoronoiSphere(DYNAMICO_reader(ncread, VoronoiMeshes(meshname)) ; prec=Float32)
pv = GeoMakie.Observable(fun.(sphere.lon_i, sphere.lat_i))

@testset "ClimFlowsPlots.jl" begin
    savefig("orthographic.png", VSPlots.plot_orthographic(sphere, pv ; colormap=:berlin))
    @test true
    savefig("2D.png", VSPlots.plot_2D(sphere, pv ; colormap=:berlin))
    @test true
    savefig("native_3D.png", VSPlots.plot_native_3D(sphere, pv ; colormap=:berlin))
    @test true
end

@testset "Spectral orthographic recording" begin
    lons, lats = 0:5:360, -90:5:90
    field = [fun(deg2rad(lon), deg2rad(lat)) for lon in lons, lat in lats]
    pv = GeoMakie.Observable(field)
    fig = SSPlots.orthographic(lons, lats, pv; colormap=:berlin, colorrange=(-1, 1))
    first_frame = copy(CairoMakie.colorbuffer(fig))

    mktempdir() do dir
        path = joinpath(dir, "spectral.mp4")
        CairoMakie.record(fig, path, (1, -1)) do sign
            pv[] = sign .* field
        end
        @test filesize(path) > 0
    end

    @test CairoMakie.colorbuffer(fig) != first_frame
end
