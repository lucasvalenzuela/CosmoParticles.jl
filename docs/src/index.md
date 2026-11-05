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
properties: vectors with one value per particle and ``3 × N`` matrices for positions and velocities. A
[`ParticleCollection`](@ref) groups several particle types. The example uses random data, so the numbers in the
outputs differ from run to run:

```jldoctest usage; filter = r"[ -]*\d+(\.\d+)?(e[+-]?\d+)?"
julia> using CosmoParticles, FillArrays

julia> dm = Particles(:dm, :id => 1:1000, :pos => 20 .* rand(3, 1000) .- 10, :vel => randn(3, 1000),
                     :mass => Fill(1.0, 1000));

julia> gas = Particles(:gas, :id => 1001:1500, :pos => 20 .* rand(3, 500) .- 10, :vel => randn(3, 500),
                       :mass => rand(500), :temp => 1e4 .* rand(500));

julia> pc = ParticleCollection(dm, gas)
ParticleCollection
dm: 1000 Particles
 id mass pos vel
gas: 500 Particles
 id mass pos temp vel
```

The rows of a matrix property can be accessed with an underscore, and functions such as [`meanpos`](@ref) compute
properties of the particles:

```jldoctest usage; filter = r"[ -]*\d+(\.\d+)?(e[+-]?\d+)?"
julia> dm.pos_x[1:3]  # the x coordinates of the first three dark matter particles
3-element Vector{Float64}:
  5.03660534992575
 -9.799353958188961
 -8.380428773084962

julia> center = meanpos(pc[:gas])  # mass-weighted mean position of the gas
3-element Vector{Float64}:
  0.4670231941202256
  0.0016508576147436433
 -0.02258312410088849
```

Particles can be selected by geometries and transformed, here keeping the particles within a distance of 5 from the
center and then moving the center to the origin:

```jldoctest usage; filter = r"[ -]*\d+(\.\d+)?(e[+-]?\d+)?"
julia> inner = filter(pc, CosmoSphere(center, 5.0))
ParticleCollection
dm: 65 Particles
 id mass pos vel
gas: 30 Particles
 id mass pos temp vel

julia> translate!(inner, -center);

julia> angmomtot(inner[:gas])  # total angular momentum of the gas in the sphere
3-element Vector{Float64}:
 6.222878954559906
 1.9490265618556184
 4.520187517576433
```

The properties of all particle types are combined lazily by `pc.all`, and [`Particles`](@ref) materializes them into
regular arrays:

```jldoctest usage; filter = r"[ -]*\d+(\.\d+)?(e[+-]?\d+)?"
julia> inner.all
all: 95 Particles
 id mass pos temp vel

julia> Particles(inner.all)
all: 95 Particles
 id mass pos temp vel
```
