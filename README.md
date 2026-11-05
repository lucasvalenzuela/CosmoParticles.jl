# CosmoParticles

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucasvalenzuela.github.io/CosmoParticles.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucasvalenzuela.github.io/CosmoParticles.jl/dev)
[![Build Status](https://github.com/lucasvalenzuela/CosmoParticles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucasvalenzuela/CosmoParticles.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lucasvalenzuela/CosmoParticles.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lucasvalenzuela/CosmoParticles.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


CosmoParticles is a Julia package that provides structs and a clean interface
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

```julia
using CosmoParticles, FillArrays

dm = Particles(:dm, :id => 1:1000, :pos => 20 .* rand(3, 1000) .- 10, :vel => randn(3, 1000),
               :mass => Fill(1.0, 1000))
gas = Particles(:gas, :id => 1001:1500, :pos => 20 .* rand(3, 500) .- 10, :vel => randn(3, 500),
                :mass => rand(500), :temp => 1e4 .* rand(500))
pc = ParticleCollection(dm, gas)

dm.pos_x                    # the x coordinates of the dark matter particles
center = meanpos(pc[:gas])  # mass-weighted mean position of the gas
inner = filter(pc, CosmoSphere(center, 5.0))  # new collection with the particles inside the sphere
translate!(inner, -center)  # move the center to the origin
angmomtot(inner[:gas])      # total angular momentum of the gas in the sphere
inner.all                   # all particle types combined lazily
Particles(inner.all)        # all particle types materialized into regular arrays
```
