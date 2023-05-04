using PolyaGammaHybridSamplers
using Documenter

DocMeta.setdocmeta!(PolyaGammaHybridSamplers, :DocTestSetup, :(using PolyaGammaHybridSamplers); recursive=true)

makedocs(;
    modules=[PolyaGammaHybridSamplers],
    authors="W.Z. Horton",
    repo="https://github.com/wzhorton/PolyaGammaHybridSamplers.jl/blob/{commit}{path}#{line}",
    sitename="PolyaGammaHybridSamplers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wzhorton.github.io/PolyaGammaHybridSamplers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wzhorton/PolyaGammaHybridSamplers.jl",
    devbranch="main",
)
