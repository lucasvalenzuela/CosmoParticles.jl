module CosmoParticles

using FastUnitfulOperations
using FillArrays
using FillArrays: AbstractFillVector
using LazyArrays
using LinearAlgebra
using Statistics
using StatsBase
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

    # operations
    delete,

    # transformations
    rotate,
    rotate!,
    translate,
    translate!,
    translate_periodic,
    translate_periodic!,
    translate_periodic_to_center,
    translate_periodic_to_center!,
    to_comoving,
    to_comoving!,
    to_physical,
    to_physical!,

    # properties
    meanprop,
    sumprop,
    meanpos,
    meanvel,
    colnorm,
    colnorm2,
    colcross,
    coldot,
    angmom,
    angmomtot,
    angmomtot_stable


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
