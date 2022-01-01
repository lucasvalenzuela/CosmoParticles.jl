module CosmoParticles

using FastUnitfulOperations
using FillArrays
using LazyArrays
using LinearAlgebra
using Unitful

export
    # types
    AbstractParticles,
    Particles,
    AllParticles,
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
    translate!,
    to_comoving,
    to_comoving!,
    to_physical,
    to_physical!,

    # properties
    colnorm,
    colnorm2,
    colcross,
    coldot


include("abstractparticles.jl")
include("particles.jl")
include("abstractparticlecollection.jl")
include("particlecollection.jl")
include("allparticles.jl")
include("geometry.jl")

include("utils.jl")
include("transformations.jl")
include("miscoperations.jl")
include("properties.jl")

end
