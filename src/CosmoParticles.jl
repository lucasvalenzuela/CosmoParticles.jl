module CosmoParticles

using FastUnitfulOperations
using LinearAlgebra
import LinearAlgebra.rotate!

export AbstractParticles, Particles, rotate, rotate!, translate, translate!

include("abstractparticles.jl")
include("particles.jl")

include("utils.jl")
include("transformations.jl")
end
