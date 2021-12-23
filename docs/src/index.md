```@meta
CurrentModule = CosmoParticles
```

# CosmoParticles

[`CosmoParticles`](https://github.com/lucasvalenzuela/CosmoParticles.jl) provides structs and a clean interface
for working with sets of particles and galaxies that have individual properties, especially made for dealing
with the data extracted from cosmological simulations.


## Installation

This package is currently only available via the local registry
[CosmoSimsRegistry](https://github.com/lucasvalenzuela/CosmoSimsRegistry),
which can be loaded from the Julia REPL with the Julia package manager:
```julia
using Pkg
pkg"registry add https://github.com/lucasvalenzuela/CosmoSimsRegistry"
```
This only needs to be done once per Julia installation.

The latest version of the package is available for Julia 1.7 and newer versions and can be installed with:
```julia
using Pkg
Pkg.add("CosmoParticles")
```
