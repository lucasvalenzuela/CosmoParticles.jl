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
translate_periodic
translate_periodic!
translate_periodic_to_center
translate_periodic_to_center!
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
delete
Base.delete!
CosmoParticles.applyind
CosmoParticles.applyind!
CosmoParticles.removeind
CosmoParticles.removeind!
CosmoParticles.findall_in
CosmoParticles.findall_in_sorted
```


## Internals

```@docs
CosmoParticles.matrix_rotate
CosmoParticles.matrix_rotate!
CosmoParticles._applyind
CosmoParticles._removeind
CosmoParticles._translate_periodic
```
