module CosmoParticles

using FastUnitfulOperations
using FillArrays
using FillArrays: AbstractFill, AbstractFillVector, AbstractFillMatrix
using LazyArrays
using LinearAlgebra
using NamedRowArrays
using NamedRowArrays: rownames
using OrderedCollections: OrderedDict
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

    # ellipsoids
    triaxiality,
    ellipticity,
    eccentricity,
    axis_lengths,

    # tables
    nrow,
    ncol

public
    # interface for particle types and collections
    get_props,
    get_particles,
    particle_name,
    particle_number,
    particle_collection,
    show_properties,

    # geometry
    geometry_enclosing_corners,
    geometry_enclosing_center,
    mask_in,
    mask_in!,
    CosmoUnionGeometry,
    CosmoIntersectGeometry,
    CosmoDiffGeometry,
    Rotated,
    Translated,
    rotation_matrix,
    rotation_matrix_inv,

    # indexing
    applyind,
    applyind!,
    removeind,
    removeind!,
    deleteat,
    findall_in,
    findall_in_sorted,

    # projections
    get_dims,
    get_losdim,

    # transformations
    matrix_rotate,
    matrix_rotate!,
    factor_to_comoving,
    factor_to_physical,

    # units
    product_preserve_type,
    product_preserve_type!,
    ustrip_lazy,
    ustrip_lazy!,
    uconvert_lazy,
    uconvert_lazy!


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
