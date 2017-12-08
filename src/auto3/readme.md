The script `get_aliases.jl` writes a file `petsc_constants_gen.c`, which is
valid C source code.  Compiling and linking it to Petsc and then running it
produces julia source code file `petsc_constants_gen.jl`.

The script reads the contents of the current `petsc_constants_gen.jl` file to
figure out what constants to get new values for.  Therefore, if you need a new
set of constants, you can get them from the `auto2/PETSc.jl` or
 `auto2/libPETSc_common.jl` file, put them into `petsc_constants_gen.jl`, then
do the steps described above to get new values for the constants.

Note that `petsc_constants_gen.jl` must have a particular structure:

```
typealias PetscTypeName JuliaType
global const CONSTNAME = value
global const CONSTNAME2 = value2
# ...
# last line must be blank
```
