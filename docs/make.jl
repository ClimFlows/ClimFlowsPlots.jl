using ClimFlowsPlots
using Documenter

DocMeta.setdocmeta!(ClimFlowsPlots, :DocTestSetup, :(using ClimFlowsPlots); recursive=true)

makedocs(;
    modules=[ClimFlowsPlots],
    authors="The ClimFlows contributors",
    sitename="ClimFlowsPlots.jl",
    format=Documenter.HTML(;
        canonical="https://ClimFlows.github.io/ClimFlowsPlots.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ClimFlows/ClimFlowsPlots.jl",
    devbranch="main",
)
