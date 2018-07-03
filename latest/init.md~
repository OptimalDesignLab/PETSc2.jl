# Initialization and Finalization

The user is responsible for initializing and finalizing Petsc.  The functions
here facilitate doing so.

The command `using PETSc2` will initialize MPI if it is not already initialized
and finalize it when Julia exits.
The user is responsible for initializing and finalizing `PETSc` itself.
Using the `Base.atexit()` function is useful for managing the finalization of
libraries.

Also note that the user is entirely responsible for finalizing `PETSc` objects
via [`PetscDestroy`](@ref) when they are no longer needed.  Unforturnately,
finalizers cannot be used because they do not guarantee order of destruction,
and they are run *after* `Base.atexit()` hooks, so it is possible `MPI` could be
finalized before the `PETSc` objects have been freed.

```@autodocs
Modules = [PETSc]
Pages = ["PETSc.jl"]
```
