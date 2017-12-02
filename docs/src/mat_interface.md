# Mat Interface

```@meta
CurrentModule = PETSc
```
This page describes an interface that is efficiently implemented for both
`Array`, `SparseMatrixCSC`, and `Petscmat`.  In some cases, `Base` functions
are extended with new methods, in other cases, new functions are defined for
the supported matrix types.  Using these functions allows some amount of
generic programming.


```@autodocs
Modules = [PETSc]
Pages = ["mat_interface.jl"]
```
