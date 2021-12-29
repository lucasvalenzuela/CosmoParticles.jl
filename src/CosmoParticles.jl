module CosmoParticles

using FastUnitfulOperations
using LinearAlgebra
using Unitful

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
    translate!,
    to_comoving,
    to_comoving!,
    to_physical,
    to_physical!


include("abstractparticles.jl")
include("particles.jl")
include("abstractparticlecollection.jl")
include("particlecollection.jl")
include("geometry.jl")

include("utils.jl")
include("transformations.jl")
include("miscoperations.jl")

end
