using GridapRobustDarcyForchheimer
using Documenter

DocMeta.setdocmeta!(GridapRobustDarcyForchheimer, :DocTestSetup, :(using GridapRobustDarcyForchheimer); recursive=true)

makedocs(;
    modules=[GridapRobustDarcyForchheimer],
    authors="Ricardo Ruiz Baier <ricardo.ruizbaier@monash.edu>",
    sitename="GridapRobustDarcyForchheimer.jl",
    format=Documenter.HTML(;
        canonical="https://ruizbaier.github.io/GridapRobustDarcyForchheimer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ruizbaier/GridapRobustDarcyForchheimer.jl",
    devbranch="main",
)
