module CosmoParticles

using FastUnitfulOperations
using LinearAlgebra
import LinearAlgebra.rotate!

export
    # types
    AbstractParticles,
    Particles,
    AbstractParticleCollection,
    ParticleCollection,
    RedshiftParticleCollection,

    # geometry
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

    # getters
    redshift,

    # transformations
    rotate,
    rotate!,
    translate,
    translate!

include("abstractparticles.jl")
include("particles.jl")
include("abstractparticlecollection.jl")
include("particlecollection.jl")
include("geometry.jl")

include("utils.jl")
include("transformations.jl")
include("miscoperations.jl")

end
