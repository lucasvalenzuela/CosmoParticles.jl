using CosmoParticles
using Documenter
using DocumenterInterLinks

DocMeta.setdocmeta!(CosmoParticles, :DocTestSetup, :(using CosmoParticles); recursive=true)

links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "NamedRowArrays" => "https://lucasvalenzuela.github.io/NamedRowArrays.jl/stable/",
    "OrderedCollections" => "https://juliacollections.github.io/OrderedCollections.jl/stable/",
    "Unitful" => "https://juliaphysics.github.io/Unitful.jl/stable/",
)

makedocs(;
    plugins=[links],
    modules=[CosmoParticles],
    checkdocs=:public,
    authors="Lucas Valenzuela",
    sitename="CosmoParticles.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucasvalenzuela.github.io/CosmoParticles.jl",
        assets=String[],
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Particles" => "particles.md",
        "Particle Collections" => "particlecollection.md",
        "Operations" => "operations.md",
        "Properties" => "properties.md",
        "Geometry" => "geometry.md",
        "Utils" => "utils.md",
    ],
)

deploydocs(; repo="github.com/lucasvalenzuela/CosmoParticles.jl", devbranch="main")
