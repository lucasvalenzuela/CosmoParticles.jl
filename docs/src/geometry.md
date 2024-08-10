```@meta
CurrentModule = CosmoParticles
```

# Geometry

To filter particles in a particular volume, such as a cube, a sphere, or a cylinder, geometries can be used.
This package provides a variety of different geometries for multiple dimensions, including the following:
- [Hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle)
  * [3D Rectangular Cuboid](https://en.wikipedia.org/wiki/Cuboid#Rectangular_cuboid)
  * [2D Rectangle](https://en.wikipedia.org/wiki/Rectangle)
- [Hypersphere](https://en.wikipedia.org/wiki/N-sphere)
  * [3D Sphere](https://en.wikipedia.org/wiki/Sphere)
  * [2D Circle](https://en.wikipedia.org/wiki/Circle)
- [Cylinder](https://en.wikipedia.org/wiki/Cylinder) (standing cylinder aligned with the z axis and arbitrary
  orientation)

```@docs
AbstractCosmoGeometry
CosmoParticles.geometry_enclosing_corners
CosmoParticles.geometry_enclosing_center
CosmoParticles.mask_in
CosmoParticles.mask_in!
CosmoHyperrectangle
CosmoCuboid
CosmoRectangle
CosmoHypercube
CosmoCube
CosmoSquare
CosmoHypersphere
CosmoHypersphere(::Integer, ::Number)
CosmoSphere
CosmoSphere(::Number)
CosmoCircle
CosmoCircle(::Number)
CosmoCylinder
CosmoCylinder(::CosmoStandingCylinder)
CosmoStandingCylinder
CosmoStandingCylinder(::CosmoCylinder)
CosmoUnionGeometry
CosmoIntersectGeometry
CosmoDiffGeometry
Base.union(::AbstractCosmoGeometry...)
Base.intersect(::AbstractCosmoGeometry...)
Base.setdiff(::Any, ::AbstractCosmoGeometry...)
Rotated
rotate(::AbstractCosmoGeometry, ::AbstractMatrix)
rotation_matrix
rotation_matrix_inv
Translated
translate(::AbstractCosmoGeometry, ::AbstractMatrix)
```
