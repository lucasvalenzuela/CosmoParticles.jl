module CosmoParticles

using FastUnitfulOperations
using LinearAlgebra
import LinearAlgebra.rotate!

export AbstractParticles,
    Particles,
    AbstractCosmoGeometry,
    CosmoHyperrectangle,
    CosmoCuboid,
    CosmoRectangle,
    CosmoHypercube,
    CosmoCube,
    CosmoSquare,
    CosmoHypersphere,
    CosmoSphere,
    CosmoCircle,
    CosmoCylinder,
    CosmoStandingCylinder,
    rotate,
    rotate!,
    translate,
    translate!

include("abstractparticles.jl")
include("particles.jl")
include("geometry.jl")

include("utils.jl")
include("transformations.jl")
include("miscoperations.jl")

end
