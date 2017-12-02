# Vec Interface

The purpose of the functions on this page is to provide a set of functions
that are implemented efficiently for both `AbstractVector` and [`PetscVec`](@ref).
In some cases `Base` functions are extended with new methods for 
[`PetscVec`](@ref), in other cases, new functions are created and defined
for both `AbstractVector` and [`PetscVec`](@ref).
Using these functions enables some amount of generic programming, particularly
when an operation is applied to either the local portion of a parallel
vector or when an operation is applied uniformly across all processors.

```@meta
CurrentModule = PETSc
```

```@autodocs
Modules = [PETSc]
Pages = ["vec_interface.jl"]
```
