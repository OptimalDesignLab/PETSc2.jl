# PETSc
This package provides thin wrappers for PETSc, as well as a few convenience functions that take advantage of multiple dispatch.

This package requires the MPI.jl package be installed.  Once it is installed you should be able to run both Julia and Petsc in parallel using MPI for all communication.  The testing verifies that PETSc can be used both serially and in parallel.

To use the package, simply put `using PETSc` at the top of your Julia source file.  The module exports the names of all the functions, as well as the Petsc data type aliases and constants such as `PETSC_DECIDE`.

In general, it is possible to run PETSC in parallel. To do so with 4 processors, do:

```
mpirun -np 4 julia ./name_of_file
```

Note that this launches 4 independent Julia processes.  They are not aware of each other using Julia's built-in parallelism, and MPI is used for all communications.  

To run in serial, do:
```
julia ./name_of_file
```

Even when running serially, the MPI.jl package must be installed.


An example of using a Krylov Sub-Space method to solve a linear system is in  `test/test_ksp.jl`, which solves a simple system with a Krylov Subspace method and compares the result with a direct solve using Julia's backslash operator.  This works in serial and in parallel.  It requires some variables declared at the top of `runtests.jl` to work.



## To do:
  * Handle Petsc error codes properly
  * Make the script for building Petsc more flexible, eg. allowing more configuration options like building BLAS or LAPCK, while ensure it remains completely autonomous (needed for Travis testing)
  * Wrap more PetscVec functions
  * Wrap more PetscMat functions
  * Wrap more KSP function
  * Determine priorities for wrapping additional functions
  * Add backup finalizer/destructor in case the variable holding the pointer to a Petsc object gets reassigned before the Petsc object is destroyed
  * Provide methods that copy indices to proper size integer
  * Create warning system (`@inefficient` macro?) for use of copying methods
  * use `Ref{T}` instead of `&` for passing by reference
  * Alias Petsc functions to Base functions where unambiguous

## Status
### Vector
Currently, most of the Vector Beginner and Intermediate functions are implimented.  Index sets are not fully supported yet.
### Matrix
 Many of the Beginner Matrix functions are implimented, although preallocation is not yet (this is coming soon)

### KSP
 Just enough KSP functions are implimented to do a GMRES solve.  After the vector and matrix functions I will focus on KSP.

## For Contributors:
  Wrappers generated by Clang.jl are in the src/auto directory.  Although not quite usable, the functions can be made useable with a few simple modifications:
  * Pass `comm.val` instead of `comm` itself for MPI communicators, and change the type to `comm_type`, which is typealiased to the the type used by the MPI.jl package
  * Pass obj.pobj instead of obj for Petsc objects as Ptr{Void}
  * For each Petsc object, you must create a type that holds a void pointer called pobj and use that in place of the (incorrect) type generated by Clang.jl
  * For every function you add, create a test

## Compatability with AbstractArray
A limited vocabulary of functions are provided that can operate on Petsc vectors
and matrices as well as AbstractArrays.  They are:

`set_values1!`
`get_values1!`
`PetscMatAssemblyBegin`
`PetscMatAssemblyEnd`


## Directory Structure
  `/src` : source files.  PETSc.jl is the main file containing initialization, with the functions for each type of Petsc object in its own file.  All constants are declared in `petsc_constants.jl`.

  `/src/auto`: auto generated wrappers from Clang.jl.  Not directly useful, but easy to modify to make useful

  `/test` : contains `runtest.jl`, which does some setup and runs all tests on the current PETSc installation.  Tests for each type of Petsc object (mirroring the files in `/src`) are contained in separate files.  The file `runtests.sh` builds PETSc and runs the tests on combinations of integer size, floating point precision, and type of scalar (real or complex).

  `/deps` : builds Petsc if needed.  See description below


## Building Petsc (or not)
Upon installation of the package, the build script (`build.jl`) checks for the environmental variables `PETSC_DIR` and `PETSC_ARCH`.  If both are present, it does nothing, otherwise is downloads and builds Petsc using the script `install_petsc.sh`.  Note that the script assumes BLAS, LAPACK and MPI (currently only MPICH is supported) are already installed.  See `.travis.yml` to see the Ubuntu packages used for testing.  When complete, it generates two files, `use_petsc.sh` and `petsc_evars`, which contains the values of `PETSC_DIR` and `PETSC_ARCH` for the Petsc installation.

  At runtime, the module checks for the presence of the environmental variables and uses them if found, otherwise it uses the values in `petsc_evars`.  This enables you to use different Petsc installations if you choose.  When the module is imported in the user code, it auto-detects the size of the Petsc integers, the precision of floats, and whether scalars are real or complex.


## Installing MPI.jl
This package requires MPI.jl, although it is not listed in the REQUIRE file because that would download the release version of MPI.jl, which does not work.  Instead, you must use the master branch.  After you have an MPI implementation installed, Pkg.build("Petsc") will install it and then Petsc, according to the description above.  If you wish to install it manually, do:
```
  Pkg.clone("MPI")
  Pkg.build("MPI")
```

Currently, only MPI implimentations where the Fortran communicator is the same as the C communictor are supported.  This is due to the current implimentation of the MPI.jl package.  MPICH is one MPI implimentation that satisfies this requirement.  Note that OPENMPI does not. 

[![Build Status](https://travis-ci.org/JaredCrean2/Petsc.svg?branch=master)](https://travis-ci.org/JaredCrean2/Petsc)
