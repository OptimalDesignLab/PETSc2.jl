# Build System

Installing `PETSc2` requires:

 * `MPI.jl`
 * `ArrayViews.jl`
 * `PETSc` itself


## Installing Other Dependencies

`MPI.jl` will be installed automatically if it is not already present.
The release version of `MPI.jl` does not work, so the build system installs
an older version. Before you install `MPI.jl`, you must have an MPI
implementation installed.  Some systems come with MPI pre-installed, others
do not.  If you are unsure which MPI implementation to install, I recommend
MPICH, which is available as a package for Debian-based system.

`ArrayViews.jl` will also be installed if not present.  


## Installing PETSc

The build system looks for an existing PETSc installation and only attempts
to build PETSc if one is not found.  It looks for the `PETSC_DIR` and
`PETSC_ARCH` environment variables.  If they are both present, PETSc is
not installed.

If you want the build system to build PETSc, you must first install

 * MPI (described in the previous section)
 * BLAS
 * LAPACK
 * Python

Note that PETSc cannot link to the BLAS and LAPACK packaged with Julia.
Instead, it links to the system BLAS and LAPACK.
The efficiency of PETSc depends on the efficiency of BLAS, so it
is recommended to obtain a high-quality BLAS (for example, if your machine has
a vendor-provided BLAS implemention) and compile it from source
on your machine.

Python is required by the PETSc build system.  It is not used by PETSc at
runtime.

By default, the build system builds a `debug` version of PETSc.  This enables
error checking within PETSc that is useful for debugging, at the cost of
some performance.  To build an optimized version of PETSc, see the `deps/install_petsc.sh` script.

Note that arguments passed to the `install_petsc.sh` script are passed to
PETSc's `configure` command.  Also, the environment variable `JULIA_PETSC_CONFIG` is added to the arguments to `configure`, before the script arguments (so
the script arguments will take precedence if an option is specified in both).
Either of these methods can be used to customize the Petsc build.

Note that you can use the `install_petsc.sh` script to install Petsc indendently
of the Petsc.jl build system.
After doing so, set the `PETSC_DIR` and `PETSC_ARCH` environment variables before
running `build.jl` (see the next section).

Even if the build system builds PETSc, it is still possible to use other
versions of PETSc by setting the `PETSC_DIR` and `PETSC_ARCH` environment
variables at runtime.  If these variables are not present, the version of
PETSc built by the build system is used.


## Linking other programs to PETSc

The build system generates a file `deps/use_petsc.sh`.  Running `source deps/use_petsc.sh` will set the `PETSC_DIR` and `PETSC_ARCH` environment variables
corresponding to the PETSc build by the build system.  This enables, for example, compiling C programs against this version of PETSc.


## Debugging Install

The PETSc library is installed in `deps/petsc-x.y.z`, where `x.y.z` is the current version of PETSc.
Inside this directory, PETSc will create a directory for your machine's specific architecture.
On my machine the architecture is `arch-linux2-c-debug`.
Inside this directory are the usual `/bin`, `/lib` etc. directories.

In the process of building PETSc, output is redirected to log files in the
`deps/petsc-x.y.z` directory.
Specifically

 * The output of `configure` is written to `fout`
 * The output of `make` is written to `fout2`
 * The output of `make test` is written to `fout3`

If any of these commands fail, the log files contain detail information that
explain what happened.
