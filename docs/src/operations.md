```@meta
CurrentModule = CosmoParticles
```

# Operations

## Transformations

```@docs
rotate
LinearAlgebra.rotate!
translate
translate!
to_comoving
to_comoving!
to_physical
to_physical!
```

## Further Operations

```@docs
Base.sort
Base.sort!
Base.filter
Base.filter!
CosmoParticles.applyind
CosmoParticles.applyind!
CosmoParticles.findall_in
CosmoParticles.findall_in_sorted
```


## Internals

```@docs
CosmoParticles.matrix_rotate
CosmoParticles.matrix_rotate!
CosmoParticles._applyind
```
