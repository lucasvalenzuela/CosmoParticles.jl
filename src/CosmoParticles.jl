module CosmoParticles

using FastUnitfulOperations
using FillArrays
using FillArrays: AbstractFill, AbstractFillVector, AbstractFillMatrix
using LazyArrays
using LinearAlgebra
using Statistics
using StatsBase
using Tables
using Tables.DataAPI: nrow, ncol
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

    # cosmo geometry
    CosmoEllipse,
    CosmoEllipsoid,
    CosmoHomoeoid2D,
    CosmoHomoeoid,

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
    colnormell,
    colnormell2,
    colcross,
    coldot,
    angmom,
    angmomtot,
    angmomtot_stable,
    
    # ellipsoids
    triaxiality,
    ellipticity,
    eccentricity,

    # tables
    nrow,
    ncol


include("abstractparticles.jl")
include("particles.jl")
include("abstractparticlecollection.jl")
include("particlecollection.jl")
include("allparticles.jl")
include("geometry.jl")
include("ellipsoids.jl")

include("utils.jl")
include("transformations.jl")
include("miscoperations.jl")
include("properties.jl")
include("tables.jl")

end
