```@meta
CurrentModule = CosmoParticles
```

# CosmoParticles

[`CosmoParticles`](https://github.com/lucasvalenzuela/CosmoParticles.jl) provides structs and a clean interface
for working with sets of particles and galaxies that have individual properties, especially made for dealing
with the data extracted from cosmological simulations.


## Installation

This package is available through the public registry
[CosmoAnalysisRegistry](https://github.com/lucasvalenzuela/CosmoAnalysisRegistry). Add the registry
once per Julia installation, together with the General registry, which provides the dependencies:
```julia
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(url="https://github.com/lucasvalenzuela/CosmoAnalysisRegistry")
```

The package requires Julia 1.12 or newer and can then be installed with:
```julia
Pkg.add("CosmoParticles")
```

## Usage

Particles of one type, such as dark matter or gas, are stored in a [`Particles`](@ref) object as a dictionary of
properties: vectors with one value per particle, ``3 × N`` matrices for positions and velocities, or a scalar for a
property that all particles share. A [`ParticleCollection`](@ref) groups several particle types:

```julia
using CosmoParticles

dm = Particles(:dm, :id => collect(1:1000), :pos => 20 .* rand(3, 1000) .- 10, :vel => randn(3, 1000), :mass => 1.0)
gas = Particles(:gas, :id => collect(1001:1500), :pos => 20 .* rand(3, 500) .- 10, :vel => randn(3, 500),
                :mass => rand(500), :temp => 1e4 .* rand(500))
pc = ParticleCollection(dm, gas)

dm.pos_x                    # the x coordinates of the dark matter particles
center = meanpos(pc[:gas])  # mass-weighted mean position of the gas
inner = filter(pc, CosmoSphere(center, 5.0))   # new collection with the particles inside the sphere
translate!(inner, -center)  # move the center to the origin
angmomtot(inner[:gas])      # total angular momentum of the gas in the sphere
inner.all.mass              # the masses of all particle types, concatenated lazily
```
