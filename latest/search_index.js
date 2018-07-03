var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "PETSc2.jl Documentation",
    "title": "PETSc2.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#PETSc2.jl-Documentation-1",
    "page": "PETSc2.jl Documentation",
    "title": "PETSc2.jl Documentation",
    "category": "section",
    "text": "This package provides an interface to the Portable, Extensible Toolkit for Scientific Computation (PETSc). Both PETSc and these wrappers are designed for solving the systems of equations that arise from discretizations of partial differential equations. Note that PETSc is not a general-purpose sparse-matrix library.Before using this any features of package, users must read the corresponding section of the PETSc manual.This package consists of two parts.  The first part is a direct wrapping of (parts of) the PETSc API.  Users should refer to the  PETSc documentation for these functions.  The documentation on this website only notes differences from the PETSc documentation.The second part is a small API that can be consistently and efficiently applied to both PETSc and regular Julia arrays (including SparseMatrixCSC)."
},

{
    "location": "index.html#First-Part-1",
    "page": "PETSc2.jl Documentation",
    "title": "First Part",
    "category": "section",
    "text": "Pages = [\"build.md\", \"vec.md\", \"mat.md\", \"constants.md\", \"ksp.md\", \"pc.md\", \"options.md\", \"error.md\"]\nDepth = 1Note that PETSc has an extensive API and not all of it has been wrapped yet. In most cases, the PETSc options database can be used rather than the API. For example, there is an extensive API for setting different options for the Krylov solver, but they can be set using the options database as follows:# Create first KSP\n\n# define options\nopts = Dict{ASCIIString, ASCIIString}(\n\"-ksp_gmres_modifiedgramschmidt => \"\"\n\"-ksp_gmres_restart\" => \"30,\n)\n\nPetscSetOptions(opts)  # set options in PETSc\'s database\nksp_one = KSP(MPI.COMM_WORLD)\nKSPSetFromOptions(ksp_one)  # copy options from PETSc\'s database to the ksp object\n  \n# Create a second KSP\n# after KSPSetFromOptions(ksp_one) is called, we are free to change the options\n# in PETSc\'s databse\nopts[\"-ksp_gmres_restart\"] = \"60\"\nPetscSetOptions(opts)  \nksp_two = KSP(MPI.COMM_WORLD)\nKSPSetFromOptions(ksp_two)The only case where PETSc\'s API is needed is when the options need to be changed as the program runs (for example, inexact Newton-Krylov requires changing the KSP solve tolerance dynamically)."
},

{
    "location": "index.html#Second-Part-1",
    "page": "PETSc2.jl Documentation",
    "title": "Second Part",
    "category": "section",
    "text": "This small API is implemented efficiently for Array, SparseMatrixCSC as well as PetscMat and PetscVec.Pages = [\"vec_interface.md\", \"mat_interface.md\"]\nDepth = 1"
},

{
    "location": "build.html#",
    "page": "Build System",
    "title": "Build System",
    "category": "page",
    "text": ""
},

{
    "location": "build.html#Build-System-1",
    "page": "Build System",
    "title": "Build System",
    "category": "section",
    "text": "Installing PETSc2 required:MPI.jl\nArrayViews.jl\nPETSc itself"
},

{
    "location": "build.html#Installing-Other-Dependencies-1",
    "page": "Build System",
    "title": "Installing Other Dependencies",
    "category": "section",
    "text": "MPI.jl will be install automatically if it is not already present. The release version of MPI.jl does not work, so the build system installs an older version. Before you install MPI.jl, you must have an MPI implementation installed.  Some systems come with MPI pre-installed, others do not.  If you are unsure which MPI implementation to install, I recommend MPICH, which is available as a package for Debian-based system.ArrayViews.jl will also be installed if not present.  "
},

{
    "location": "build.html#Installing-PETSc-1",
    "page": "Build System",
    "title": "Installing PETSc",
    "category": "section",
    "text": "The build system looks for an existing PETSc installation and only attempts to build PETSc if one is not found.  It looks for the PETSC_DIR and PETSC_ARCH environment variables.  If they are both present, PETSc is not installed.If you want the build system to build PETSc, you must first installMPI (described in the previous section)\nBLAS\nLAPACK\nPythonNote that PETSc cannot link to the BLAS and LAPACK packaged with Julia. Instead, it links to the system BLAS and LAPACK. The efficiency of PETSc depends on the efficiency of BLAS, so it is recommended to obtain a high-quality BLAS (for example, if you machine has a vendor-provided BLAS implemention) and compile it from source on your machine.Python is required by the PETSc build system.  It is not used by PETSc at runtime.By default, the build system builds a debug version of PETSc.  This enables error checking within PETSc that is useful for debugging, at the cost of some performance.  To build an optimized version of PETSc, see the deps/install_petsc.sh script.Note that arguments passed to the install_petsc.sh script are passed to PETSc\'s configure command.  If you want non-default options, running this script directly (not as part of the build system) is the recommended way to build PETSc. After doing so, set the PETSC_DIR and PETSC_ARCH environment variables before running build.jl (see the next section).Even if the build system builds PETSc, it is still possible to use other versions of PETSc by setting the PETSC_DIR and PETSC_ARCH environment variables at runtime.  If these variables are not present, the version of PETSc built by the build system is used."
},

{
    "location": "build.html#Linking-other-programs-to-PETSc-1",
    "page": "Build System",
    "title": "Linking other programs to PETSc",
    "category": "section",
    "text": "The build system generates a file deps/use_petsc.sh.  Running source deps/use_petsc.sh will set the PETSC_DIR and PETSC_ARCH environment variables corresponding to the PETSc build by the build system.  This enables, for example, compiling C programs against this version of PETSc."
},

{
    "location": "build.html#Debugging-Install-1",
    "page": "Build System",
    "title": "Debugging Install",
    "category": "section",
    "text": "The PETSc library is installed in deps/petsc-x.y.z, where x.y.z is the current version of PETSc. Inside this directory, PETSc will create a directory for your machine\'s specific architecture. On my machine the architecture is arch-linux2-c-debug. Inside this directory are the usual /bin, /lib etc. directories.In the process of building PETSc, output is redirected to log files in the deps/petsc-x.y.z directory. SpecificallyThe output of configure is written to fout\nThe output of make is written to fout2\nThe output of make test is written to fout3If any of these commands fail, the log files contain detail information that explain what happened."
},

{
    "location": "init.html#",
    "page": "Initialization and Finalization",
    "title": "Initialization and Finalization",
    "category": "page",
    "text": ""
},

{
    "location": "init.html#Initialization-and-Finalization-1",
    "page": "Initialization and Finalization",
    "title": "Initialization and Finalization",
    "category": "section",
    "text": "The user is responsible for initializing and finalizing Petsc.  The functions here facilitate doing so.The command using PETSc2 will initialize MPI if it is not already initialized and finalize it when Julia exits. The user is responsible for initializing and finalizing PETSc itself. Using the Base.atexit() function is useful for managing the finalization of libraries.Also note that the user is entirely responsible for finalizing PETSc objects via PetscDestroy when they are no longer needed.  Unforturnately, finalizers cannot be used because they do not guarantee order of destruction, and they are run after Base.atexit() hooks, so it is possible MPI could be finalized before the PETSc objects have been freed.Modules = [PETSc2]\nPages = [\"PETSc.jl\"]"
},

{
    "location": "vec.html#",
    "page": "Vec Documentation",
    "title": "Vec Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "vec.html#PETSc2.PetscVec",
    "page": "Vec Documentation",
    "title": "PETSc2.PetscVec",
    "category": "type",
    "text": "Petsc Vector type.\n\nNot a subtype of AbstractArray because Petsc vectors do not conform to that API.\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.PetscDestroy-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "Free a Petsc vec.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.PetscView-Tuple{PETSc2.PetscVec,Any}",
    "page": "Vec Documentation",
    "title": "PETSc2.PetscView",
    "category": "method",
    "text": "PetscView for Petsc vector\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAXPBY-Tuple{PETSc2.PetscVec,Float64,Float64,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAXPBY",
    "category": "method",
    "text": "VecAXPBY\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAXPBYPCZ-Tuple{PETSc2.PetscVec,Float64,Float64,Float64,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAXPBYPCZ",
    "category": "method",
    "text": "VecAXPBYPCZ\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAXPY-Tuple{PETSc2.PetscVec,Float64,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAXPY",
    "category": "method",
    "text": "VecAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAYPX-Tuple{PETSc2.PetscVec,Float64,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAYPX",
    "category": "method",
    "text": "VecAYPX\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAssemblyBegin-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAssemblyBegin",
    "category": "method",
    "text": "VecAssemblyBegin\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAssemblyEnd-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAssemblyEnd",
    "category": "method",
    "text": "VecAssemblyEnd\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecCopy-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecCopy",
    "category": "method",
    "text": "VecCopy\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecDot-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecDot",
    "category": "method",
    "text": "VecDot\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecDuplicate-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecDuplicate",
    "category": "method",
    "text": "VecDulicate\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecExp-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecExp",
    "category": "method",
    "text": "VecExp\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetArray-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetArray",
    "category": "method",
    "text": "VecGetArray.  Users must call VecRestoreArray when finished.\n\nInputs\n\nvec: the Petsc vector\n\nOutputs\n\narr: a Julia Array{PetscScalar, 1}\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetArrayRead-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetArrayRead",
    "category": "method",
    "text": "Similar to VecGetArray, but produces the array must not be written to.\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetLocalSize-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetLocalSize",
    "category": "method",
    "text": "VecGetLocalSize\n\nInputs\n\nvec: a Petsc vector\n\nOutputs\n\nthe size\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetOwnershipRange-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetOwnershipRange",
    "category": "method",
    "text": "VecGetOwnershipRange\n\nInputs\n\nvec: Petsc vector\n\nOutputs\n\nlow: lowest index (zero-based) that is owned\nhigh: highest index + 1 that is owned\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetSize-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetSize",
    "category": "method",
    "text": "VeGetSize\n\nInputs\n\nvec: a Petsc vector\n\nOutputs\n\nthe size\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetValues-Tuple{PETSc2.PetscVec,AbstractArray{Int32,1},AbstractArray{Float64,1}}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetValues",
    "category": "method",
    "text": "VecGetValues with length of idx inferred from idx\n\nInputs\n\nvec: the Petsc vector\nidx: array of PetscInt indices\ny: array of PetscScalar values\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecGetValues-Tuple{PETSc2.PetscVec,Integer,AbstractArray{Int32,1},AbstractArray{Float64,1}}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecGetValues",
    "category": "method",
    "text": "VecGetValues\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecLog-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecLog",
    "category": "method",
    "text": "VecLog\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecMAXPY-Tuple{PETSc2.PetscVec,Integer,AbstractArray{Float64,1},AbstractArray{Ptr{Void},1}}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecMAXPY",
    "category": "method",
    "text": "VecMAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecMax-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecMax",
    "category": "method",
    "text": "VecMax\n\nInputs\n\nvec: PetscVec\n\nOutput\n\nr: the maximum value\nidx: the (zero-based) index of the maximum value\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecMin-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecMin",
    "category": "method",
    "text": "VecMin. Same interface as VecMax\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecNorm-Tuple{PETSc2.PetscVec,Integer}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecNorm",
    "category": "method",
    "text": "VecNorm\n\nInputs\n\nobj: Petsc vector\nnormtype: the Petsc enum for the norm type\n\nOutput\n\nthe norm value\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecPointwiseDivide-Tuple{PETSc2.PetscVec,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecPointwiseDivide",
    "category": "method",
    "text": "VecPointwiseDivide\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecPointwiseMult-Tuple{PETSc2.PetscVec,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecPointwiseMult",
    "category": "method",
    "text": "VecPointwiseMult\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecReciprocal-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecReciprocal",
    "category": "method",
    "text": "VecReciprocal\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecRestoreArray-Tuple{PETSc2.PetscVec,Array{Float64,1}}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecRestoreArray",
    "category": "method",
    "text": "VecRestoreArray.  Users must not access the array after calling this function.\n\nInputs\n\nvec: the PetscVector passed into VecGetArray\narr: the array returned by `VecGetArray\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecRestoreArrayRead-Tuple{PETSc2.PetscVec,Array{Float64,1}}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecRestoreArrayRead",
    "category": "method",
    "text": "Similar to VecRestoreArray, but corresponds to VecGetArrayRead\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecScale-Tuple{PETSc2.PetscVec,Float64}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecScale",
    "category": "method",
    "text": "VecScale\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSet-Tuple{PETSc2.PetscVec,Float64}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSet",
    "category": "method",
    "text": "VecSet\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSetSizes-Tuple{PETSc2.PetscVec,Integer,Integer}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSetSizes",
    "category": "method",
    "text": "VecSetSizes\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSetType-Tuple{PETSc2.PetscVec,Any}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSetType",
    "category": "method",
    "text": "VecSetType\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSetValues-Tuple{PETSc2.PetscVec,Array{Int32,N} where N,Array{Float64,N} where N,Integer}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSetValues",
    "category": "method",
    "text": "VecSetValues\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSetValues-Tuple{PETSc2.PetscVec,Array{Int32,N} where N,Array{Float64,N} where N}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSetValues",
    "category": "method",
    "text": "VecSetValues method that implicitly uses INSERT_VALUES\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecShift-Tuple{PETSc2.PetscVec,Float64}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecShift",
    "category": "method",
    "text": "VecShift\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSqrtAbs-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSqrtAbs",
    "category": "method",
    "text": "VecSqrtAbs\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSum-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSum",
    "category": "method",
    "text": "VecSum\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecSwap-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecSwap",
    "category": "method",
    "text": "VecSwap\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecTDot-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecTDot",
    "category": "method",
    "text": "VecTDot\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecWAXPY-Tuple{PETSc2.PetscVec,Float64,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecWAXPY",
    "category": "method",
    "text": "VecWAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.AllVectors",
    "page": "Vec Documentation",
    "title": "PETSc2.AllVectors",
    "category": "constant",
    "text": "Union{AbstractVector, PetscVec}\n\n\n\n"
},

{
    "location": "vec.html#PETSc2.VecAssemble-Tuple{PETSc2.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc2.VecAssemble",
    "category": "method",
    "text": "Convenience function for calling VecAssemblyBegin, and VecAssemblyEnd, in   one go.\n\n\n\n"
},

{
    "location": "vec.html#Vec-Documentation-1",
    "page": "Vec Documentation",
    "title": "Vec Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page wraps functions in PETSc\'s Vec API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc2]\nPages = [\"vec.jl\"]"
},

{
    "location": "vec_interface.html#",
    "page": "Vec Interface",
    "title": "Vec Interface",
    "category": "page",
    "text": ""
},

{
    "location": "vec_interface.html#PETSc2.PetscVec-Tuple{Integer,Any,MPI.Comm}",
    "page": "Vec Interface",
    "title": "PETSc2.PetscVec",
    "category": "method",
    "text": "Create a PetscVec, setting both the type and the format.  Users can specify   either the local or global dimensions\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.PetscVec-Tuple{Integer,MPI.Comm}",
    "page": "Vec Interface",
    "title": "PETSc2.PetscVec",
    "category": "method",
    "text": "Create a vector of a given size.  Users can specify either the global   dimension or the local dimension\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.PetscDestroy-Tuple{AbstractArray{T,1} where T}",
    "page": "Vec Interface",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "PetscDestroy for AbstractVector.  No-op\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.PetscView",
    "page": "Vec Interface",
    "title": "PETSc2.PetscView",
    "category": "function",
    "text": "Print non-Petsc vector to a given IO (a Julia IO, not a Petsc IO).  Defaults   to printing to STDOUT.\n\nInputs\n\nb: AbstractVector\nf: IO, defaults to STDOUT\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.assembly_begin-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.assembly_begin",
    "category": "method",
    "text": "Calls VecAssemblyBegin.  No-op for Julia vectors\n\nInputs\n\nvec: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.assembly_end-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.assembly_end",
    "category": "method",
    "text": "Calls VecAssemblyEnd.  No-op for Julia vectors\n\nInputs\n\nvec: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.diagonal_shift!-Tuple{PETSc2.PetscVec,Number}",
    "page": "Vec Interface",
    "title": "PETSc2.diagonal_shift!",
    "category": "method",
    "text": "Add a given value to all elements of the vector\n\nInputs\n\nA: AbstractVector\na: number to shift by\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.fill_zero!-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.fill_zero!",
    "category": "method",
    "text": "Fill a vector with zeros\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.get_values1!-Tuple{AbstractArray{T,1} where T,Array{Int32,N} where N,Array}",
    "page": "Vec Interface",
    "title": "PETSc2.get_values1!",
    "category": "method",
    "text": "Method for AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.get_values1!-Tuple{PETSc2.PetscVec,Array{Int32,N} where N,Array{Float64,N} where N}",
    "page": "Vec Interface",
    "title": "PETSc2.get_values1!",
    "category": "method",
    "text": "Like set_values1! but for retrieving values.  Note that Petsc   only supports retrieving values from the local part of the vector\n\nInputs\n\nvec: a vector, can be a julia vector or a Petsc vector.\n\nInputs/Outputs\n\nidx: indices to retrieve\nvals: array to put the values into\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.length_global-Tuple{Union{AbstractArray{T,1} where T, PETSc2.PetscVec}}",
    "page": "Vec Interface",
    "title": "PETSc2.length_global",
    "category": "method",
    "text": "Length of global vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.length_local-Tuple{Union{AbstractArray{T,1} where T, PETSc2.PetscVec}}",
    "page": "Vec Interface",
    "title": "PETSc2.length_local",
    "category": "method",
    "text": "Length of local part of vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.local_indices-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.local_indices",
    "category": "method",
    "text": "Returns a UnitRange containing the (1-based) global indices owned by this   process.\n\nInputs\n\nA: AbstractVector\n\nOutputs\n\nrng: UnitRange\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.set_values1!",
    "page": "Vec Interface",
    "title": "PETSc2.set_values1!",
    "category": "function",
    "text": "1-based indexing for both regular vectors and Petsc vector\n\nInputs\n\nvals: values to add/insert into the vector, must be length(idx)\nflag: INSERT_VALUES or ADD_VALUES\n\nInputs/Outputs\n\nvec: the vector, can be a Petsc vector or a julia vector\nidx: (global) indices to add/insert vals into\n\nidx is listed as input/output because it may be modified during the function.   It will be returned to its original values when the function exits.   This is necessary to accomodate Petscs zero-based indexing interface\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.set_values1!-Union{Tuple{AbstractArray{T,1} where T,Array{Int32,N} where N,Array{T,N} where N,Integer}, Tuple{AbstractArray{T,1} where T,Array{Int32,N} where N,Array{T,N} where N}, Tuple{T}} where T",
    "page": "Vec Interface",
    "title": "PETSc2.set_values1!",
    "category": "method",
    "text": "Method for AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.size_global-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.size_global",
    "category": "method",
    "text": "Size of global vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc2.size_local-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc2.size_local",
    "category": "method",
    "text": "Size of local part of vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.dot-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.LinAlg.dot",
    "category": "method",
    "text": "Base.dot\n\nDot product where the first vector is conjugated.  This is is the reverse   of VecDot, where the second vector is conjugated\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.norm",
    "page": "Vec Interface",
    "title": "Base.LinAlg.norm",
    "category": "function",
    "text": "Base.norm\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.scale!-Tuple{PETSc2.PetscVec,Number}",
    "page": "Vec Interface",
    "title": "Base.LinAlg.scale!",
    "category": "method",
    "text": "Base.scale!\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.copy!-Tuple{PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.copy!",
    "category": "method",
    "text": "Base.copy!\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.copy-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.copy",
    "category": "method",
    "text": "Base.copy\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.fill!-Tuple{PETSc2.PetscVec,Any}",
    "page": "Vec Interface",
    "title": "Base.fill!",
    "category": "method",
    "text": "Base.fill! for Petsc vectors\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.maximum-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.maximum",
    "category": "method",
    "text": "Base.maximum\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.minimum-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.minimum",
    "category": "method",
    "text": "Base.minimum\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.sum-Tuple{PETSc2.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.sum",
    "category": "method",
    "text": "Base.sum\n\n\n\n"
},

{
    "location": "vec_interface.html#Vec-Interface-1",
    "page": "Vec Interface",
    "title": "Vec Interface",
    "category": "section",
    "text": "The purpose of the functions on this page is to provide a set of functions that are implemented efficiently for both AbstractVector and PetscVec. In some cases Base functions are extended with new methods for  PetscVec, in other cases new functions are created and defined for both AbstractVector and PetscVec. Using these functions enables some amount of generic programming, particularly when an operation is applied to either the local portion of a parallel vector or when an operation is applied uniformly across all processors.CurrentModule = PETSc2Modules = [PETSc2]\nPages = [\"vec_interface.jl\"]"
},

{
    "location": "mat.html#",
    "page": "Mat Documentation",
    "title": "Mat Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "mat.html#PETSc2.PetscMat",
    "page": "Mat Documentation",
    "title": "PETSc2.PetscMat",
    "category": "type",
    "text": "PetscMat type.  Currently a subtype of AbstractArray, although that may   change.\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAXPY-Tuple{PETSc2.PetscMat,Float64,PETSc2.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAXPY",
    "category": "method",
    "text": "MatAxPY\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAYPX-Tuple{PETSc2.PetscMat,Float64,PETSc2.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAYPX",
    "category": "method",
    "text": "MatAYPX\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAssemblyBegin-Tuple{PETSc2.PetscMat,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAssemblyBegin",
    "category": "method",
    "text": "MatAssemblyBegin\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAssemblyBegin-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAssemblyBegin",
    "category": "method",
    "text": "MatAssemblyBegin, impicitly using MAT_FINAL_ASSEMBLY\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAssemblyEnd-Tuple{PETSc2.PetscMat,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAssemblyEnd",
    "category": "method",
    "text": "MatAssemblyEnd\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatAssemblyEnd-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatAssemblyEnd",
    "category": "method",
    "text": "MatAssemblyEnd, implicitly using MAT_FINAL_ASSEMBLY\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatCreateShell-Tuple{MPI.Comm,Integer,Integer,Integer,Integer,Ptr{Void}}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatCreateShell",
    "category": "method",
    "text": "MatCreateShell\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatCreateTranspose-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatCreateTranspose",
    "category": "method",
    "text": "MatCreateTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetInfo-Tuple{PETSc2.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetInfo",
    "category": "method",
    "text": "MatGetInfo\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetLocalSize-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetLocalSize",
    "category": "method",
    "text": "MatGetLocalsize.  Same interface as MatGetSize\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetOwnershipRange-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetOwnershipRange",
    "category": "method",
    "text": "MatGetOwnershipRange\n\nInputs\n\nmat: a PetscMat\n\nOutputs\n\nlow: lowest (zero-based) index that is owned\nhigh: highest + 1 (zero-based) index that is owned\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetSize-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetSize",
    "category": "method",
    "text": "MatGetSize\n\nInputs\n\nobj: a PetscMat\n\nOutputs\n\nm: first dimension size\nn: second dimension size\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetType-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetType",
    "category": "method",
    "text": "MatGetType\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatGetValues-Tuple{PETSc2.PetscMat,Array{Int32,1},Array{Int32,1},Array{Float64,2}}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatGetValues",
    "category": "method",
    "text": "MatGetValues\n\nInputs\n\nobj: the PetscMat\nidxm: row indices\nidxn: column indices\nvals: logically 2 dimensional array of values (although it could be a       vector too)\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMPIAIJSetPreallocation-Tuple{PETSc2.PetscMat,Integer,Union{AbstractArray{Int32,N} where N, Ptr{Void}},Integer,Union{AbstractArray{Int32,N} where N, Ptr{Void}}}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMPIAIJSetPreallocation",
    "category": "method",
    "text": "MatMPIAIJSetPreallocation\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMatMult-Tuple{PETSc2.PetscMat,PETSc2.PetscMat,Int32,Float64,PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMatMult",
    "category": "method",
    "text": "MatMatMult\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMult-Tuple{PETSc2.PetscMat,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMult",
    "category": "method",
    "text": "MatMult\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMultAdd-Tuple{PETSc2.PetscMat,PETSc2.PetscVec,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMultAdd",
    "category": "method",
    "text": "MatMultAdd\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMultHermitianTranspose-Tuple{PETSc2.PetscMat,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMultHermitianTranspose",
    "category": "method",
    "text": "MatMultHermitianTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatMultTranspose-Tuple{PETSc2.PetscMat,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatMultTranspose",
    "category": "method",
    "text": "MatMultTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatNorm-Tuple{PETSc2.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatNorm",
    "category": "method",
    "text": "MatNorm\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatScale-Tuple{PETSc2.PetscMat,Number}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatScale",
    "category": "method",
    "text": "MatScale\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatSetOption-Tuple{PETSc2.PetscMat,Int32,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatSetOption",
    "category": "method",
    "text": "MatSetOption\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatSetSizes-Tuple{PETSc2.PetscMat,Integer,Integer,Integer,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatSetSizes",
    "category": "method",
    "text": "MatSetSizes\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatSetType-Tuple{PETSc2.PetscMat,Any}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatSetType",
    "category": "method",
    "text": "MatSetType\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatSetValues-Tuple{PETSc2.PetscMat,Array{Int32,N} where N,Array{Int32,N} where N,Array{Float64,N} where N,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatSetValues",
    "category": "method",
    "text": "MatSetValues.  The idx and idy arrays must have PetscInt elements, and the   vals array must have PetscScalar elements.\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatSetValuesBlocked-Tuple{PETSc2.PetscMat,Array{Int32,N} where N,Array{Int32,N} where N,Array{Float64,N} where N,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatSetValuesBlocked",
    "category": "method",
    "text": "MatSetValues blocked.\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatShellGetContext-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatShellGetContext",
    "category": "method",
    "text": "MatShellGetContext\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatShellSetOperation-Tuple{PETSc2.PetscMat,UInt32,Ptr{Void}}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatShellSetOperation",
    "category": "method",
    "text": "MatShellSetOperation\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatShift-Tuple{PETSc2.PetscMat,Number}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatShift",
    "category": "method",
    "text": "MatShift\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatTranspose-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatTranspose",
    "category": "method",
    "text": "Constructs the transpose of a given matrix.  Can be in place or out of place.\n\nInputs\n\nA: matrix to take the transpose of\n\nKeywords\n\ninplace: whether or not to do the transpose in place\n\nOutputs\n\nA PetscMat object, either the original object A if the transpose was done  in place, or a new matrix object if the transpose was done out of place\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatXAIJSetPreallocation-Tuple{PETSc2.PetscMat,Integer,AbstractArray{Int32,1},AbstractArray{Int32,1},AbstractArray{Int32,1},AbstractArray{Int32,1}}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatXAIJSetPreallocation",
    "category": "method",
    "text": "MatXAIJSetPreallocation\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatZeroEntries-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.MatZeroEntries",
    "category": "method",
    "text": "MatZeroEntries\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.PetscDestroy-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "Frees a Petsc object.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.PetscView",
    "page": "Mat Documentation",
    "title": "PETSc2.PetscView",
    "category": "function",
    "text": "PetscViewer\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.SetFromOptions-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.SetFromOptions",
    "category": "method",
    "text": "MatSetFromOptions\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.SetUp-Tuple{PETSc2.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc2.SetUp",
    "category": "method",
    "text": "MatSetUp\n\n\n\n"
},

{
    "location": "mat.html#PETSc2.MatInfo",
    "page": "Mat Documentation",
    "title": "PETSc2.MatInfo",
    "category": "type",
    "text": "Equivalent to Petsc\'s MatInfo struct\n\n\n\n"
},

{
    "location": "mat.html#Base.show-Tuple{IO,PETSc2.MatInfo}",
    "page": "Mat Documentation",
    "title": "Base.show",
    "category": "method",
    "text": "show() for a MatInfo\n\n\n\n"
},

{
    "location": "mat.html#Mat-Documentation-1",
    "page": "Mat Documentation",
    "title": "Mat Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page wraps functions in PETSc\'s Mat API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc2]\nPages = [\"mat.jl\"]"
},

{
    "location": "mat_interface.html#",
    "page": "Mat Interface",
    "title": "Mat Interface",
    "category": "page",
    "text": ""
},

{
    "location": "mat_interface.html#PETSc2.PetscMat-Tuple{Integer,Integer,Any,MPI.Comm}",
    "page": "Mat Interface",
    "title": "PETSc2.PetscMat",
    "category": "method",
    "text": "Constructor\n\nInputs\n\nmglobal: first dimension global size (or PETSC_DECIDE)\nnglobal: second dimension global size (or PETSC_DECIDE)\nformat: matrix format\ncomm: MPI communicator\n\nKeyword Arguments\n\nmlocal: first dimension local size (or PETSC_DECIDE)\nnlocal: second dimension local size (or PETSC_DECIDE)\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.PetscDestroy-Tuple{AbstractArray{T,2} where T}",
    "page": "Mat Interface",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "PetscDestroy for AbstractMatrix.  No-op\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.PetscView",
    "page": "Mat Interface",
    "title": "PETSc2.PetscView",
    "category": "function",
    "text": "Print a non-Petsc matrix to a given IO (a Julia IO, not a Petsc IO)\n\nInputs\n\nA: the matrix\nf: an IO, default STDOUT\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.assembly_begin-Tuple{PETSc2.PetscMat,Integer}",
    "page": "Mat Interface",
    "title": "PETSc2.assembly_begin",
    "category": "method",
    "text": "Begin matrix assembly for PetscMat.  No-op for Julia matrices\n\nInputs\n\nmat: AbstractMatrix\nflag: type of matrix assembly (see MatAssemblyBegin\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.assembly_end-Tuple{PETSc2.PetscMat,Integer}",
    "page": "Mat Interface",
    "title": "PETSc2.assembly_end",
    "category": "method",
    "text": "Counterpart of assembly_end\n\nInputs\n\nmat: AbstractMatrix\nflg: type of matrix assembly\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.diagonal_shift!-Tuple{PETSc2.PetscMat,Number}",
    "page": "Mat Interface",
    "title": "PETSc2.diagonal_shift!",
    "category": "method",
    "text": "Adds the specified value to the diagonal of the matrix\n\nInputs\n\nA: the matrix\na: the value\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.fill_zero!-Tuple{PETSc2.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc2.fill_zero!",
    "category": "method",
    "text": "Fill the matrix with zeros. The sparsity pattern of the matrix (if applicable)   should be defined before this function is called\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.get_values1!-Tuple{AbstractArray{T,2} where T,Array{Int32,N} where N,Array{Int32,N} where N,Array}",
    "page": "Mat Interface",
    "title": "PETSc2.get_values1!",
    "category": "method",
    "text": "Method for AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.get_values1!-Tuple{PETSc2.PetscMat,Array{Int32,N} where N,Array{Int32,N} where N,Array{Float64,N} where N}",
    "page": "Mat Interface",
    "title": "PETSc2.get_values1!",
    "category": "method",
    "text": "Like set_values1!, but retrieves values.  See that function for   the meanings of the arguments. Note that Petsc does   not support getting values for the non-local block of the matrix\n\nInputs\n\nmat: a matrix, can be a Petsc matrix or a julia matrix\n\nInputs/Outputs\n\nidxm\nidxn\nvals\n\nAliasing restrictions: idxm and idxn cannot alias\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.local_indices-Tuple{PETSc2.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc2.local_indices",
    "category": "method",
    "text": "Returns the rows owned by this process (1-based)\n\nInputs\n\nA: a matrix\n\nOutputs\n\nrng: a UnitRange\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.set_values1!",
    "page": "Mat Interface",
    "title": "PETSc2.set_values1!",
    "category": "function",
    "text": "1-based indexing for both regular and Pets matrices.   Note that Petsc treats arrays at being row-major, so it is recommened   to set MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE) before using   this function.\n\nInputs\n\nflag: INSERT_VALUES or ADD_VALUES.  Note that the first one       result in non-deterministic behavior in parallel (in the general       case)\nvals: the values, must be length(idxm) x length(idxn)\n\nInputs/Outputs\n\nmat: a matrix, can be a Petsc matrix or a julia matrix\nidxm: the row numbers\nidxn: the column numbers\n\nNote that idxm and idxn are listed as input/outputs because they may be   modified by this function, however when the function returns they   will have the same values as on entry.  This is needed to accomodate the   fact that Petsc uses 1 based indexing internally.\n\nThis function is optimized for PetscMat and SparseMatrixCSC\n\nAliasing restriction: idxm and idxn cannot alias\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.set_values1!-Union{Tuple{AbstractArray{T,2} where T,Array{Int32,N} where N,Array{Int32,N} where N,Array{T,N} where N,Integer}, Tuple{AbstractArray{T,2} where T,Array{Int32,N} where N,Array{Int32,N} where N,Array{T,N} where N}, Tuple{T}} where T",
    "page": "Mat Interface",
    "title": "PETSc2.set_values1!",
    "category": "method",
    "text": "Method for AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.set_values1!-Union{Tuple{SparseMatrixCSC,Array{Int32,N} where N,Array{Int32,N} where N,Array{T,N} where N,Integer}, Tuple{SparseMatrixCSC,Array{Int32,N} where N,Array{Int32,N} where N,Array{T,N} where N}, Tuple{T}} where T",
    "page": "Mat Interface",
    "title": "PETSc2.set_values1!",
    "category": "method",
    "text": "Method for SparseMatrixCSC\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.size_global-Tuple{PETSc2.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc2.size_global",
    "category": "method",
    "text": "Global size of matrix, same as size_local() for serial matrices\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc2.size_local-Tuple{PETSc2.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc2.size_local",
    "category": "method",
    "text": "Size of local part of matrix\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.A_mul_B!-Tuple{PETSc2.PetscVec,PETSc2.PetscMat,PETSc2.PetscVec}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.A_mul_B!",
    "category": "method",
    "text": "Computes b = A*x\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.At_mul_B!-Tuple{PETSc2.PetscVec,PETSc2.PetscMat,PETSc2.PetscVec}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.At_mul_B!",
    "category": "method",
    "text": "Computes b = A.\'*x\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.norm-Tuple{PETSc2.PetscMat,Real}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.norm",
    "category": "method",
    "text": "Norm for Petsc matrices, 1, 2, and infinity norms supported\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.scale!-Tuple{PETSc2.PetscMat,Number}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.scale!",
    "category": "method",
    "text": "scale! for Petsc matrix\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.vecnorm-Tuple{PETSc2.PetscMat}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.vecnorm",
    "category": "method",
    "text": "Frobenius norm, consistent with Julias interface\n\n\n\n"
},

{
    "location": "mat_interface.html#Mat-Interface-1",
    "page": "Mat Interface",
    "title": "Mat Interface",
    "category": "section",
    "text": "CurrentModule = PETSc2This page describes an interface that is efficiently implemented for both Array, SparseMatrixCSC, and Petscmat.  In some cases, Base functions are extended with new methods, in other cases new functions are defined for the supported matrix types.  Using these functions allows some amount of generic programming.Modules = [PETSc2]\nPages = [\"mat_interface.jl\"]"
},

{
    "location": "constants.html#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "constants.html#PETSc2.PetscInt",
    "page": "Constants",
    "title": "PETSc2.PetscInt",
    "category": "type",
    "text": "Petsc integer type\n\n\n\n"
},

{
    "location": "constants.html#PETSc2.PetscReal",
    "page": "Constants",
    "title": "PETSc2.PetscReal",
    "category": "type",
    "text": "The closest real type to PetscScalar\n\n\n\n"
},

{
    "location": "constants.html#PETSc2.PetscScalar",
    "page": "Constants",
    "title": "PETSc2.PetscScalar",
    "category": "type",
    "text": "Element type of Petsc vectors and matrices\n\n\n\n"
},

{
    "location": "constants.html#Constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": "CurrentModule = PETSc2Many PETSc constants and enums are available in the Julia wrappers.  The full list can be found in src/petsc_constants.jl.  See the PETSc documentation for their meaning. A few of the most important constants are listed here.  Not all of the constants are exported from the PETSc module, to avoid cluttering the user\'s namespace.Modules = [PETSc2]\nPages = [\"constants.jl\"]"
},

{
    "location": "ksp.html#",
    "page": "KSP Documentation",
    "title": "KSP Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "ksp.html#PETSc2.KSP_NULL",
    "page": "KSP Documentation",
    "title": "PETSc2.KSP_NULL",
    "category": "constant",
    "text": "Null pointer KSP object\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.KSP",
    "page": "KSP Documentation",
    "title": "PETSc2.KSP",
    "category": "type",
    "text": "KSP object\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetConvergedReason-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetConvergedReason",
    "category": "method",
    "text": "KSPGetConvergedReason\n\nInputs\n\nksp: KSP object\n\nOutputs\n\nstring containing reason\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetInitialGuessNonzero-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetInitialGuessNonzero",
    "category": "method",
    "text": "KSPGetInitialGuessNonzero\n\nInputs\n\nKSP\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetOperators-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetOperators",
    "category": "method",
    "text": "KSPGetOperators\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetResidualNorm-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetResidualNorm",
    "category": "method",
    "text": "KSPGetResidualNorm\n\nInputs\n\nKSP\n\nOutputs\n\nPetscReal\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetTolerances-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetTolerances",
    "category": "method",
    "text": "KSPGetTolerances\n\nInputs\n\nKSP\n\nOutputs\n\nrtol\nabstol\ndtol\nmaxits\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.GetType-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.GetType",
    "category": "method",
    "text": "KSPGetType\n\nInputs\n\nKSP\n\nOutputs\n\nstring containing KSP type\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.KSPSolve-Tuple{PETSc2.KSP,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "KSP Documentation",
    "title": "PETSc2.KSPSolve",
    "category": "method",
    "text": "KSPSolve\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.KSPSolveTranspose-Tuple{PETSc2.KSP,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "KSP Documentation",
    "title": "PETSc2.KSPSolveTranspose",
    "category": "method",
    "text": "KSPSolveTranspose\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.PetscDestroy-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "PetscDestroy for KSP object.  Safe to call multiple times.\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.PetscView",
    "page": "KSP Documentation",
    "title": "PETSc2.PetscView",
    "category": "function",
    "text": "PetscView\n\nInputs\n\nksp: KSP object\nviewer: PetscViewer, defaults to Petsc stdout\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetFromOptions-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetFromOptions",
    "category": "method",
    "text": "SetFromOptions\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetInitialGuessNonzero-Tuple{PETSc2.KSP,UInt32}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetInitialGuessNonzero",
    "category": "method",
    "text": "KSPSetInitialGuessNonzero\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetOperators-Tuple{PETSc2.KSP,PETSc2.PetscMat,PETSc2.PetscMat}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetOperators",
    "category": "method",
    "text": "KSPSetOperators\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetTolerances-Tuple{PETSc2.KSP,Number,Number,Number,Integer}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetTolerances",
    "category": "method",
    "text": "KSPSetTolerances\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetType-Tuple{PETSc2.KSP,String}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetType",
    "category": "method",
    "text": "KSPSetType\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetUp-Tuple{PETSc2.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetUp",
    "category": "method",
    "text": "KSPSetUp\n\n\n\n"
},

{
    "location": "ksp.html#PETSc2.SetReusePreconditioner-Tuple{PETSc2.KSP,UInt32}",
    "page": "KSP Documentation",
    "title": "PETSc2.SetReusePreconditioner",
    "category": "method",
    "text": "KSPSetReusePreconditioner\n\n\n\n"
},

{
    "location": "ksp.html#KSP-Documentation-1",
    "page": "KSP Documentation",
    "title": "KSP Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page wraps functions in PETSc\'s Vec API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc2]\nPages = [\"ksp.jl\"]"
},

{
    "location": "pc.html#",
    "page": "PC Documentation",
    "title": "PC Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "pc.html#PETSc2.PC",
    "page": "PC Documentation",
    "title": "PETSc2.PC",
    "category": "type",
    "text": "Petsc PC object\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PC-Tuple{MPI.Comm}",
    "page": "PC Documentation",
    "title": "PETSc2.PC",
    "category": "method",
    "text": "Constructor\n\nInputs\n\ncomm: MPI communicator\n\nOutputs\n\nPC object\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.KSPGetPC-Tuple{PETSc2.KSP}",
    "page": "PC Documentation",
    "title": "PETSc2.KSPGetPC",
    "category": "method",
    "text": "KSPGetPC\n\nInputs\n\nksp: KSP object\n\nOutput\n\nPC: pc object\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.KSPSetPC-Tuple{PETSc2.KSP,PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.KSPSetPC",
    "category": "method",
    "text": "KSPSetPC\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCApply-Tuple{PETSc2.PC,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "PC Documentation",
    "title": "PETSc2.PCApply",
    "category": "method",
    "text": "PCApply\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCApplyTranspose-Tuple{PETSc2.PC,PETSc2.PetscVec,PETSc2.PetscVec}",
    "page": "PC Documentation",
    "title": "PETSc2.PCApplyTranspose",
    "category": "method",
    "text": "PCApplyTranspose\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCApplyTransposeExists-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCApplyTransposeExists",
    "category": "method",
    "text": "PCApplyTransposeExists\n\nInputs\n\npc: PC object\n\nOutputs\n\nBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCBJacobiGetSubKSP-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCBJacobiGetSubKSP",
    "category": "method",
    "text": "PCBJacobiGetSubKSP\n\nInputs\n\npc: PC object\n\nOutputs\n\nn_local: number of local KSP object\nfirst_local: global number of first KSP object on this block\nksp_arr: array of KSP object\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorGetAllowDiagonalFill-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorGetAllowDiagonalFill",
    "category": "method",
    "text": "PCFactorGetAllowDiagonalFill\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorGetLevels-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorGetLevels",
    "category": "method",
    "text": "PCFactorGetLevels\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscInt\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorGetUseInPlace-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorGetUseInPlace",
    "category": "method",
    "text": "PCFactorGetUseInPlace\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorSetAllowDiagonalFill-Tuple{PETSc2.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorSetAllowDiagonalFill",
    "category": "method",
    "text": "PCFactorSetAllowDiagonalFill\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorSetFill-Tuple{PETSc2.PC,Float64}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorSetFill",
    "category": "method",
    "text": "PCFactorSetFill\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorSetLevels-Tuple{PETSc2.PC,Int32}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorSetLevels",
    "category": "method",
    "text": "PCFactorSetLevels\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCFactorSetUseInPlace-Tuple{PETSc2.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc2.PCFactorSetUseInPlace",
    "category": "method",
    "text": "PCFactorSetUseInPlace\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCGetReusePreconditioner-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCGetReusePreconditioner",
    "category": "method",
    "text": "PCGetReusePreconditioner\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCGetType-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCGetType",
    "category": "method",
    "text": "PCGetType\n\nInputs\n\npc: PC object\n\nOutputs\n\nstring contianing the PC type\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCJacobiGetType-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCJacobiGetType",
    "category": "method",
    "text": "PCJacobiGetType\n\nInputs\n\npc: PC object\n\nOutputs\n\nstring containing the type\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCJacobiSetType-Tuple{PETSc2.PC,Int32}",
    "page": "PC Documentation",
    "title": "PETSc2.PCJacobiSetType",
    "category": "method",
    "text": "PCJacobiaSetType\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCSetReusePreconditioner-Tuple{PETSc2.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc2.PCSetReusePreconditioner",
    "category": "method",
    "text": "PCSetReusePreconditioner\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCSetType-Tuple{PETSc2.PC,String}",
    "page": "PC Documentation",
    "title": "PETSc2.PCSetType",
    "category": "method",
    "text": "PCSetType\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCSetUp-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCSetUp",
    "category": "method",
    "text": "PCSetUp\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCShellGetContext-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PCShellGetContext",
    "category": "method",
    "text": "PCShellGetContext\n\nInputs\n\npc: PC object\n\nOutputs\n\nPtr{Void}.  Users should call unsafe_pointer_to_objref() on it to get the             Julia object back\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCShellSetApply-Tuple{PETSc2.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc2.PCShellSetApply",
    "category": "method",
    "text": "PCShellSetApply\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCShellSetApplyTranspose-Tuple{PETSc2.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc2.PCShellSetApplyTranspose",
    "category": "method",
    "text": "PCShellSetApplyTranspose\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCShellSetContext-Tuple{PETSc2.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc2.PCShellSetContext",
    "category": "method",
    "text": "PCShellSetContext\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PCShellSetSetUp-Tuple{PETSc2.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc2.PCShellSetSetUp",
    "category": "method",
    "text": "PCShellSetSetUp\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.PetscDestroy-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.PetscDestroy",
    "category": "method",
    "text": "Free a Petsc PC object.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "pc.html#PETSc2.SetFromOptions-Tuple{PETSc2.PC}",
    "page": "PC Documentation",
    "title": "PETSc2.SetFromOptions",
    "category": "method",
    "text": "PCSetFromOptions\n\n\n\n"
},

{
    "location": "pc.html#PC-Documentation-1",
    "page": "PC Documentation",
    "title": "PC Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page wraps functions in PETSc\'s PC API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc2]\nPages = [\"pc.jl\"]"
},

{
    "location": "options.html#",
    "page": "Options Documentation",
    "title": "Options Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "options.html#PETSc2.PetscClearOptions-Tuple{Dict}",
    "page": "Options Documentation",
    "title": "PETSc2.PetscClearOptions",
    "category": "method",
    "text": "Convenience wrapper for using a dictionary to clear options (only the keys   are used).\n\n\n\n"
},

{
    "location": "options.html#PETSc2.PetscOptionsClearValue",
    "page": "Options Documentation",
    "title": "PETSc2.PetscOptionsClearValue",
    "category": "function",
    "text": "PetscOptionsClearValue\n\nInputs\n\narg1: the key (string)\narg2: the PetscOptions object, defaults to the global options databse\n\n\n\n"
},

{
    "location": "options.html#PETSc2.PetscOptionsSetValue",
    "page": "Options Documentation",
    "title": "PETSc2.PetscOptionsSetValue",
    "category": "function",
    "text": "PetscOptionsSetValue\n\nInputs\n\narg1: the key (string)\narg2: the value (string)\narg3: the PetscOptions object, defaults to the global options database\n\n\n\n"
},

{
    "location": "options.html#PETSc2.PetscOptionsView",
    "page": "Options Documentation",
    "title": "PETSc2.PetscOptionsView",
    "category": "function",
    "text": "PetscOptionsView\n\nInputs\n\narg1: a PetscViewer, defaults to Petsc\'s stdout\narg2: the PetscOptions object, defaults to the global options databse\n\n\n\n"
},

{
    "location": "options.html#PETSc2.PetscSetOptions-Tuple{Dict}",
    "page": "Options Documentation",
    "title": "PETSc2.PetscSetOptions",
    "category": "method",
    "text": "Convenience wrapper for using a dictionary to set options\n\n\n\n"
},

{
    "location": "options.html#PETSc2.PetscOptions",
    "page": "Options Documentation",
    "title": "PETSc2.PetscOptions",
    "category": "type",
    "text": "Typedef of PetscOptions\n\n\n\n"
},

{
    "location": "options.html#Options-Documentation-1",
    "page": "Options Documentation",
    "title": "Options Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page wraps functions in PETSc\'s options API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc2]\nPages = [\"options.jl\"]"
},

{
    "location": "error.html#",
    "page": "Error Handling Documentation",
    "title": "Error Handling Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "error.html#PETSc2.error_handler-Tuple{MPI.CComm,Int32,Ptr{UInt8},Ptr{UInt8},Int32,Int32,Ptr{UInt8},Ptr{Void}}",
    "page": "Error Handling Documentation",
    "title": "PETSc2.error_handler",
    "category": "method",
    "text": "Error handler registered with PETSc.  This function prints the information   Petsc supplies about the error, a Julia stack trace, and then invokes the   PetscTraceBackErrorHandler to print the Petsc stack trace.\n\n\n\n"
},

{
    "location": "error.html#Error-Handling-Documentation-1",
    "page": "Error Handling Documentation",
    "title": "Error Handling Documentation",
    "category": "section",
    "text": "CurrentModule = PETSc2This page describes error handling in the Julia wrappers.Modules = [PETSc2]\nPages = [\"error.jl\"]"
},

]}
