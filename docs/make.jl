using CosmoParticles
using Documenter

DocMeta.setdocmeta!(CosmoParticles, :DocTestSetup, :(using CosmoParticles); recursive=true)

makedocs(;
    modules=[CosmoParticles],
    authors="Lucas Valenzuela",
    repo="https://github.com/lucasvalenzuela/CosmoParticles.jl/blob/{commit}{path}#{line}",
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

deploydocs(;
    repo="github.com/lucasvalenzuela/CosmoParticles.jl",
    devbranch="main",
)
