# PETSc.jl Documentation

This package provides an interface to the Portable, Extensible Toolkit for
Scientific Computation ([PETSc](https://www.mcs.anl.gov/petsc/)).
Both PETSc and these wrappers are designed for solving the systems of equations
that arise from discretizations of partial differential equations.
Note that PETSc is not a general-purpose sparse-matrix library.

Before using this any features of package, users *must* read the corresponding
section of the [PETSc manual](http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf).

This package consists of two parts.  The first part is a direct wrapping of
(parts of) the PETSc API.  Users should refer to the 
[PETSc documentation](https://www.mcs.anl.gov/petsc/documentation/index.html)
for these functions.  The documentation on this website only notes differences
from the PETSc documentation.

The second part is a small API that can be consistently and efficiently
applied to both PETSc and regular Julia arrays (including `SparseMatrixCSC`).

## First Part

```@contents
Pages = ["vec.md", "mat.md", "constants.md", "ksp.md", "pc.md", "options.md", "error.md"]
Depth = 1
```

Note that PETSc has an extensive API and not all of it has been wrapped yet.
In most cases, the PETSc options database can be used rather than the API.
For example, there is an [extensive API](http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/index.html) for setting different options for the
Krylov solver, but they can be set using the options database as follows:

```
# Create first KSP

# define options
opts = Dict{ASCIIString, ASCIIString}(
"-ksp_gmres_modifiedgramschmidt => ""
"-ksp_gmres_restart" => "30,
)

PetscSetOptions(opts)  # set options in PETSc's database
ksp_one = KSP(MPI.COMM_WORLD)
KSPSetFromOptions(ksp_one)  # copy options from PETSc's database to the ksp object
  
# Create a second KSP
# after KSPSetFromOptions(ksp_one) is called, we are free to change the options
# in PETSc's databse
opts["-ksp_gmres_restart"] = "60"
PetscSetOptions(opts)  
ksp_two = KSP(MPI.COMM_WORLD)
KSPSetFromOptions(ksp_two)
```

The only case where PETSc's API is needed is when the options need to be changed
as the program runs (for example, inexact Newton-Krylov requires changing the
KSP solve tolerance dynamically).

## Second Part

This small API is implemented efficiently for `Array`, `SparseMatrixCSC` as well
as `PetscMat` and `PetscVec`.

```@contents
Pages = ["vec_interface.md", "mat_interface.md"]
Depth = 1
```
