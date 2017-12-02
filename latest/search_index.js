var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "PETSc.jl Documentation",
    "title": "PETSc.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#PETSc.jl-Documentation-1",
    "page": "PETSc.jl Documentation",
    "title": "PETSc.jl Documentation",
    "category": "section",
    "text": "This package provides an interface to the Portable, Extensible Toolkit for Scientific Computation (PETSc). Both PETSc and these wrappers are designed for solving the systems of equations that arise from discretizations of partial differential equations. Note that PETSc is not a general-purpose sparse-matrix library.Before using this any features of package, users must read the corresponding section of the PETSc manual.This package consists of two parts.  The first part is a direct wrapping of (parts of) the PETSc API.  Users should refer to the  PETSc documentation for these functions.  The documentation on this website only notes differences from the PETSc documentation.The second part is a small API that can be consistently and efficiently applied to both PETSc and regular Julia arrays (including SparseMatrixCSC)."
},

{
    "location": "index.html#First-Part-1",
    "page": "PETSc.jl Documentation",
    "title": "First Part",
    "category": "section",
    "text": "Pages = [\"vec.md\", \"mat.md\", \"constants.md\", \"ksp.md\", \"pc.md\", \"options.md\", \"error.md\"]\nDepth = 1Note that PETSc has an extensive API and not all of it has been wrapped yet. In most cases, the PETSc options database can be used rather than the API. For example, there is an extensive API for setting different options for the Krylov solver, but they can be set using the options database as follows:# Create first KSP\n\n# define options\nopts = Dict{ASCIIString, ASCIIString}(\n\"-ksp_gmres_modifiedgramschmidt => \"\"\n\"-ksp_gmres_restart\" => \"30,\n)\n\nPetscSetOptions(opts)  # set options in PETSc's database\nksp_one = KSP(MPI.COMM_WORLD)\nKSPSetFromOptions(ksp_one)  # copy options from PETSc's database to the ksp object\n  \n# Create a second KSP\n# after KSPSetFromOptions(ksp_one) is called, we are free to change the options\n# in PETSc's databse\nopts[\"-ksp_gmres_restart\"] = \"60\"\nPetscSetOptions(opts)  \nksp_two = KSP(MPI.COMM_WORLD)\nKSPSetFromOptions(ksp_two)The only case where PETSc's API is needed is when the options need to be changed as the program runs (for example, inexact Newton-Krylov requires changing the KSP solve tolerance dynamically)."
},

{
    "location": "index.html#Second-Part-1",
    "page": "PETSc.jl Documentation",
    "title": "Second Part",
    "category": "section",
    "text": "This small API is implemented efficiently for Array, SparseMatrixCSC as well as PetscMat and PetscVec.Pages = [\"vec_interface.md\", \"mat_interface.md\"]\nDepth = 1"
},

{
    "location": "init.html#",
    "page": "Initialization and Finalization",
    "title": "Initialization and Finalization",
    "category": "page",
    "text": ""
},

{
    "location": "init.html#PETSc.PetscFinalize-Tuple{}",
    "page": "Initialization and Finalization",
    "title": "PETSc.PetscFinalize",
    "category": "Method",
    "text": "PetscFinalize\n\n\n\n"
},

{
    "location": "init.html#PETSc.PetscInitialize-Tuple{Any}",
    "page": "Initialization and Finalization",
    "title": "PETSc.PetscInitialize",
    "category": "Method",
    "text": "PetscInitialize\n\nInputs\n\nargs: Array{ASCIIString, 1} containing Petsc options keys and values\n\n\n\n"
},

{
    "location": "init.html#PETSc.PetscInitialize-Tuple{}",
    "page": "Initialization and Finalization",
    "title": "PETSc.PetscInitialize",
    "category": "Method",
    "text": "PetscInitialize\n\n\n\n"
},

{
    "location": "init.html#PETSc.PetscInitialized-Tuple{}",
    "page": "Initialization and Finalization",
    "title": "PETSc.PetscInitialized",
    "category": "Method",
    "text": "PetscInitialized\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "init.html#Initialization-and-Finalization-1",
    "page": "Initialization and Finalization",
    "title": "Initialization and Finalization",
    "category": "section",
    "text": "The user is responsible for initializing and finalizing Petsc.  The functions here facilitate doing so. It is recommended for users to initialize MPI before executing using PETSc, and then finalize MPI when finished (after finalizing PETSc). Otherwise, PETSc will initialize MPI and finalize it when Julia exits. Using the Base.atexit() function is useful for managing the finalization of libraries.Also note that the user is entirely responsible for finalizing PETSc objects via PetscDestroy when they are no longer needed.  Unforturnately, finalizers cannot be used because they do not guarantee order of destruction and they are run after Base.atexit() hooks, so it is possible MPI could be finalized before the PETSc objects have been freed.Modules = [PETSc]\nPages = [\"PETSc.jl\"]"
},

{
    "location": "vec.html#",
    "page": "Vec Documentation",
    "title": "Vec Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "vec.html#PETSc.PetscVec",
    "page": "Vec Documentation",
    "title": "PETSc.PetscVec",
    "category": "Type",
    "text": "Petsc Vector type.\n\nNot a subtype of AbstractArray because Petsc vectors do not conform to that API.\n\n\n\n"
},

{
    "location": "vec.html#PETSc.PetscDestroy-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "Free a Petsc vec.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "vec.html#PETSc.PetscView-Tuple{PETSc.PetscVec,Any}",
    "page": "Vec Documentation",
    "title": "PETSc.PetscView",
    "category": "Method",
    "text": "PetscView for Petsc vector\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAXPBY-Tuple{PETSc.PetscVec,Complex{Float64},Complex{Float64},PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAXPBY",
    "category": "Method",
    "text": "VecAXPBY\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAXPBYPCZ-Tuple{PETSc.PetscVec,Complex{Float64},Complex{Float64},Complex{Float64},PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAXPBYPCZ",
    "category": "Method",
    "text": "VecAXPBYPCZ\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAXPY-Tuple{PETSc.PetscVec,Complex{Float64},PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAXPY",
    "category": "Method",
    "text": "VecAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAYPX-Tuple{PETSc.PetscVec,Complex{Float64},PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAYPX",
    "category": "Method",
    "text": "VecAYPX\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAssemblyBegin-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAssemblyBegin",
    "category": "Method",
    "text": "VecAssemblyBegin\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAssemblyEnd-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAssemblyEnd",
    "category": "Method",
    "text": "VecAssemblyEnd\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecCopy-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecCopy",
    "category": "Method",
    "text": "VecCopy\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecDot-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecDot",
    "category": "Method",
    "text": "VecDot\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecDuplicate-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecDuplicate",
    "category": "Method",
    "text": "VecDulicate\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecExp-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecExp",
    "category": "Method",
    "text": "VecExp\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetArray-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetArray",
    "category": "Method",
    "text": "VecGetArray.  Users must call VecRestoreArray when finished.\n\nInputs\n\nvec: the Petsc vector\n\nOutputs\n\narr: a Julia Array{PetscScalar, 1}\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetArrayRead-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetArrayRead",
    "category": "Method",
    "text": "Similar to VecGetArray, but produces the array must not be written to.\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetLocalSize-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetLocalSize",
    "category": "Method",
    "text": "VecGetLocalSize\n\nInputs\n\nvec: a Petsc vector\n\nOutputs\n\nthe size\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetOwnershipRange-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetOwnershipRange",
    "category": "Method",
    "text": "VecGetOwnershipRange\n\nInputs\n\nvec: Petsc vector\n\nOutputs\n\nlow: lowest index (zero-based) that is owned\nhigh: highest index + 1 that is owned\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetSize-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetSize",
    "category": "Method",
    "text": "VeGetSize\n\nInputs\n\nvec: a Petsc vector\n\nOutputs\n\nthe size\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetValues-Tuple{PETSc.PetscVec,AbstractArray{Int64,1},AbstractArray{Complex{Float64},1}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetValues",
    "category": "Method",
    "text": "VecGetValues with length of idx inferred from idx\n\nInputs\n\nvec: the Petsc vector\nidx: array of PetscInt indices\ny: array of PetscScalar values\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecGetValues-Tuple{PETSc.PetscVec,Integer,AbstractArray{Int64,1},AbstractArray{Complex{Float64},1}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecGetValues",
    "category": "Method",
    "text": "VecGetValues\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecLog-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecLog",
    "category": "Method",
    "text": "VecLog\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecMAXPY-Tuple{PETSc.PetscVec,Integer,AbstractArray{Complex{Float64},1},AbstractArray{Ptr{Void},1}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecMAXPY",
    "category": "Method",
    "text": "VecMAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecMax-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecMax",
    "category": "Method",
    "text": "VecMax\n\nInputs\n\nvec: PetscVec\n\nOutput\n\nr: the maximum value\nidx: the (zero-based) index of the maximum value\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecMin-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecMin",
    "category": "Method",
    "text": "VecMin. Same interface as VecMax\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecNorm-Tuple{PETSc.PetscVec,Integer}",
    "page": "Vec Documentation",
    "title": "PETSc.VecNorm",
    "category": "Method",
    "text": "VecNorm\n\nInputs\n\nobj: Petsc vector\nnormtype: the Petsc enum for the norm type\n\nOutput\n\nthe norm value\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecPointwiseDivide-Tuple{PETSc.PetscVec,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecPointwiseDivide",
    "category": "Method",
    "text": "VecPointwiseDivide\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecPointwiseMult-Tuple{PETSc.PetscVec,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecPointwiseMult",
    "category": "Method",
    "text": "VecPointwiseMult\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecReciprocal-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecReciprocal",
    "category": "Method",
    "text": "VecReciprocal\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecRestoreArray-Tuple{PETSc.PetscVec,Array{Complex{Float64},1}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecRestoreArray",
    "category": "Method",
    "text": "VecRestoreArray.  Users must not access the array after calling this function.\n\nInputs\n\nvec: the PetscVector passed into \nVecGetArray\narr: the array returned by `VecGetArray\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecRestoreArrayRead-Tuple{PETSc.PetscVec,Array{Complex{Float64},1}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecRestoreArrayRead",
    "category": "Method",
    "text": "Similar to VecRestoreArray, but corresponds to VecGetArrayRead\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecScale-Tuple{PETSc.PetscVec,Complex{Float64}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecScale",
    "category": "Method",
    "text": "VecScale\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSet-Tuple{PETSc.PetscVec,Complex{Float64}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSet",
    "category": "Method",
    "text": "VecSet\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSetSizes-Tuple{PETSc.PetscVec,Integer,Integer}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSetSizes",
    "category": "Method",
    "text": "VecSetSizes\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSetType-Tuple{PETSc.PetscVec,Any}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSetType",
    "category": "Method",
    "text": "VecSetType\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSetValues-Tuple{PETSc.PetscVec,Array{Int64,N},Array{Complex{Float64},N},Integer}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSetValues",
    "category": "Method",
    "text": "VecSetValues\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSetValues-Tuple{PETSc.PetscVec,Array{Int64,N},Array{Complex{Float64},N}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSetValues",
    "category": "Method",
    "text": "VecSetValues method that implicitly uses PETSC_INSERT_VALUES\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecShift-Tuple{PETSc.PetscVec,Complex{Float64}}",
    "page": "Vec Documentation",
    "title": "PETSc.VecShift",
    "category": "Method",
    "text": "VecShift\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSqrtAbs-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSqrtAbs",
    "category": "Method",
    "text": "VecSqrtAbs\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSum-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSum",
    "category": "Method",
    "text": "VecSum\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecSwap-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecSwap",
    "category": "Method",
    "text": "VecSwap\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecTDot-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecTDot",
    "category": "Method",
    "text": "VecTDot\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecWAXPY-Tuple{PETSc.PetscVec,Complex{Float64},PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecWAXPY",
    "category": "Method",
    "text": "VecWAXPY\n\n\n\n"
},

{
    "location": "vec.html#PETSc.AllVectors",
    "page": "Vec Documentation",
    "title": "PETSc.AllVectors",
    "category": "Constant",
    "text": "Union{AbstractVector, PetscVec}\n\n\n\n"
},

{
    "location": "vec.html#PETSc.VecAssemble-Tuple{PETSc.PetscVec}",
    "page": "Vec Documentation",
    "title": "PETSc.VecAssemble",
    "category": "Method",
    "text": "Convenience function for calling VecAssemblyBegin, and VecAssemblyEnd, in   one go.\n\n\n\n"
},

{
    "location": "vec.html#Vec-Documentation-1",
    "page": "Vec Documentation",
    "title": "Vec Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page wraps functions in PETSc's Vec API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc]\nPages = [\"vec.jl\"]"
},

{
    "location": "vec_interface.html#",
    "page": "Vec Interface",
    "title": "Vec Interface",
    "category": "page",
    "text": ""
},

{
    "location": "vec_interface.html#PETSc.PetscVec-Tuple{Integer,Any,MPI.Comm}",
    "page": "Vec Interface",
    "title": "PETSc.PetscVec",
    "category": "Method",
    "text": "Create a PetscVec, setting both the type and the format.  Users can specify   either the local or global dimensions\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.PetscVec-Tuple{Integer,MPI.Comm}",
    "page": "Vec Interface",
    "title": "PETSc.PetscVec",
    "category": "Method",
    "text": "Create a vector of a given size.  Users can specify either the global   dimension or the local dimension\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.PetscDestroy-Tuple{AbstractArray{T,1}}",
    "page": "Vec Interface",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "PetscDestroy for AbstractVector.  No-op\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.PetscView",
    "page": "Vec Interface",
    "title": "PETSc.PetscView",
    "category": "Function",
    "text": "Print non-Petsc vector to a given IO (a Julia IO, not a Petsc IO).  Defaults   to printing to STDOUT.\n\nInputs\n\nb: AbstractVector\nf: IO, defaults to STDOUT\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.assembly_begin-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.assembly_begin",
    "category": "Method",
    "text": "Calls VecAssemblyBegin.  No-op for Julia vectors\n\nInputs\n\nvec: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.assembly_end-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.assembly_end",
    "category": "Method",
    "text": "Calls VecAssemblyEnd.  No-op for Julia vectors\n\nInputs\n\nvec: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.diagonal_shift!-Tuple{PETSc.PetscVec,Number}",
    "page": "Vec Interface",
    "title": "PETSc.diagonal_shift!",
    "category": "Method",
    "text": "Add a given value to all elements of the vector\n\nInputs\n\nA: AbstractVector\na: number to shift by\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.fill_zero!-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.fill_zero!",
    "category": "Method",
    "text": "Fill a vector with zeros\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.get_values1!-Tuple{AbstractArray{T,1},Array{Int64,N},Array{T,N}}",
    "page": "Vec Interface",
    "title": "PETSc.get_values1!",
    "category": "Method",
    "text": "Method for AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.get_values1!-Tuple{PETSc.PetscVec,Array{Int64,N},Array{Complex{Float64},N}}",
    "page": "Vec Interface",
    "title": "PETSc.get_values1!",
    "category": "Method",
    "text": "Like set_values1! but for retrieving values.  Note that Petsc   only supports retrieving values from the local part of the vector\n\nInputs\n\nvec: a vector, can be a julia vector or a Petsc vector.   \n\nInputs/Outputs\n\nidx: indices to retrieve\nvals: array to put the values into\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.length_global-Tuple{Union{AbstractArray{T,1},PETSc.PetscVec}}",
    "page": "Vec Interface",
    "title": "PETSc.length_global",
    "category": "Method",
    "text": "Length of global vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.length_local-Tuple{Union{AbstractArray{T,1},PETSc.PetscVec}}",
    "page": "Vec Interface",
    "title": "PETSc.length_local",
    "category": "Method",
    "text": "Length of local part of vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.local_indices-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.local_indices",
    "category": "Method",
    "text": "Returns a UnitRange containing the (1-based) global indices owned by this   process.\n\nInputs\n\nA: AbstractVector\n\nOutputs\n\nrng: UnitRange\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.set_values1!",
    "page": "Vec Interface",
    "title": "PETSc.set_values1!",
    "category": "Function",
    "text": "Method for AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.set_values1!",
    "page": "Vec Interface",
    "title": "PETSc.set_values1!",
    "category": "Function",
    "text": "1-based indexing for both regular vectors and Petsc vector\n\nInputs\n\nvals: values to add/insert into the vector, must be length(idx)\nflag: PETSC_INSERT_VALUES or PETSC_ADD_VALUES\n\nInputs/Outputs\n\nvec: the vector, can be a Petsc vector or a julia vector\nidx: (global) indices to add/insert vals into      idx is listed as input/output because it may be modified during the function.   It will be returned to its original values when the function exits.   This is necessary to accomodate Petscs zero-based indexing interface\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.size_global-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.size_global",
    "category": "Method",
    "text": "Size of global vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#PETSc.size_local-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "PETSc.size_local",
    "category": "Method",
    "text": "Size of local part of vector\n\nInputs\n\nA: AbstractVector\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.dot-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.LinAlg.dot",
    "category": "Method",
    "text": "Base.dot\n\nDot product where the first vector is conjugated.  This is is the reverse   of VecDot, where the second vector is conjugated\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.norm",
    "page": "Vec Interface",
    "title": "Base.LinAlg.norm",
    "category": "Function",
    "text": "Base.norm\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.LinAlg.scale!-Tuple{PETSc.PetscVec,Number}",
    "page": "Vec Interface",
    "title": "Base.LinAlg.scale!",
    "category": "Method",
    "text": "Base.scale!\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.copy!-Tuple{PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.copy!",
    "category": "Method",
    "text": "Base.copy!\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.copy-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.copy",
    "category": "Method",
    "text": "Base.copy\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.fill!-Tuple{PETSc.PetscVec,Any}",
    "page": "Vec Interface",
    "title": "Base.fill!",
    "category": "Method",
    "text": "Base.fill! for Petsc vectors\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.maximum-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.maximum",
    "category": "Method",
    "text": "Base.maximum\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.minimum-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.minimum",
    "category": "Method",
    "text": "Base.minimum\n\n\n\n"
},

{
    "location": "vec_interface.html#Base.sum-Tuple{PETSc.PetscVec}",
    "page": "Vec Interface",
    "title": "Base.sum",
    "category": "Method",
    "text": "Base.sum\n\n\n\n"
},

{
    "location": "vec_interface.html#Vec-Interface-1",
    "page": "Vec Interface",
    "title": "Vec Interface",
    "category": "section",
    "text": "The purpose of the functions on this page is to provide a set of functions that are implemented efficiently for both AbstractVector and PetscVec. In some cases Base functions are extended with new methods for  PetscVec, in other cases, new functions are created and defined for both AbstractVector and PetscVec. Using these functions enables some amount of generic programming, particularly when an operation is applied to either the local portion of a parallel vector or when an operation is applied uniformly across all processors.CurrentModule = PETScModules = [PETSc]\nPages = [\"vec_interface.jl\"]"
},

{
    "location": "mat.html#",
    "page": "Mat Documentation",
    "title": "Mat Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "mat.html#PETSc.PetscMat",
    "page": "Mat Documentation",
    "title": "PETSc.PetscMat",
    "category": "Type",
    "text": "PetscMat type.  Currently a subtype of AbstractArray, although that may   change.\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAXPY-Tuple{PETSc.PetscMat,Complex{Float64},PETSc.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAXPY",
    "category": "Method",
    "text": "MatAxPY\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAYPX-Tuple{PETSc.PetscMat,Complex{Float64},PETSc.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAYPX",
    "category": "Method",
    "text": "MatAYPX\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAssemblyBegin-Tuple{PETSc.PetscMat,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAssemblyBegin",
    "category": "Method",
    "text": "MatAssemblyBegin\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAssemblyBegin-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAssemblyBegin",
    "category": "Method",
    "text": "MatAssemblyBegin, impicitly using PETSC_MAT_FINAL_ASSEMBLY\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAssemblyEnd-Tuple{PETSc.PetscMat,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAssemblyEnd",
    "category": "Method",
    "text": "MatAssemblyEnd\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatAssemblyEnd-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatAssemblyEnd",
    "category": "Method",
    "text": "MatAssemblyEnd, implicitly using PETSC_MAT_FINAL_ASSEMBLY\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatCreateShell-Tuple{MPI.Comm,Integer,Integer,Integer,Integer,Ptr{Void}}",
    "page": "Mat Documentation",
    "title": "PETSc.MatCreateShell",
    "category": "Method",
    "text": "MatCreateShell\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatCreateTranspose-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatCreateTranspose",
    "category": "Method",
    "text": "MatCreateTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetInfo-Tuple{PETSc.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetInfo",
    "category": "Method",
    "text": "MatGetInfo\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetLocalSize-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetLocalSize",
    "category": "Method",
    "text": "MatGetLocalsize.  Same interface as MatGetSize\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetOwnershipRange-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetOwnershipRange",
    "category": "Method",
    "text": "MatGetOwnershipRange\n\nInputs\n\nmat: a PetscMat\n\nOutputs\n\nlow: lowest (zero-based) index that is owned\nhigh: highest + 1 (zero-based) index that is owned\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetSize-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetSize",
    "category": "Method",
    "text": "MatGetSize\n\nInputs\n\nobj: a PetscMat\n\nOutputs\n\nm: first dimension size\nn: second dimension size\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetType-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetType",
    "category": "Method",
    "text": "MatGetType\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatGetValues-Tuple{PETSc.PetscMat,Array{Int64,1},Array{Int64,1},Array{Complex{Float64},2}}",
    "page": "Mat Documentation",
    "title": "PETSc.MatGetValues",
    "category": "Method",
    "text": "MatGetValues\n\nInputs\n\nobj: the PetscMat\nidxm: row indices\nidxn: column indices\nvals: logically 2 dimensional array of values (although it could be a            vector too)\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMPIAIJSetPreallocation-Tuple{PETSc.PetscMat,Integer,Union{AbstractArray{Int64,N},Ptr{Void}},Integer,Union{AbstractArray{Int64,N},Ptr{Void}}}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMPIAIJSetPreallocation",
    "category": "Method",
    "text": "MatMPIAIJSetPreallocation\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMatMult-Tuple{PETSc.PetscMat,PETSc.PetscMat,Int32,Float64,PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMatMult",
    "category": "Method",
    "text": "MatMatMult\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMult-Tuple{PETSc.PetscMat,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMult",
    "category": "Method",
    "text": "MatMult\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMultAdd-Tuple{PETSc.PetscMat,PETSc.PetscVec,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMultAdd",
    "category": "Method",
    "text": "MatMultAdd\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMultHermitianTranspose-Tuple{PETSc.PetscMat,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMultHermitianTranspose",
    "category": "Method",
    "text": "MatMultHermitianTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatMultTranspose-Tuple{PETSc.PetscMat,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "Mat Documentation",
    "title": "PETSc.MatMultTranspose",
    "category": "Method",
    "text": "MatMultTranspose\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatNorm-Tuple{PETSc.PetscMat,Int32}",
    "page": "Mat Documentation",
    "title": "PETSc.MatNorm",
    "category": "Method",
    "text": "MatNorm\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatScale-Tuple{PETSc.PetscMat,Number}",
    "page": "Mat Documentation",
    "title": "PETSc.MatScale",
    "category": "Method",
    "text": "MatScale\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetFromOptions-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetFromOptions",
    "category": "Method",
    "text": "MatSetFromOptions\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetOption-Tuple{PETSc.PetscMat,Int32,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetOption",
    "category": "Method",
    "text": "MatSetOption\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetSizes-Tuple{PETSc.PetscMat,Integer,Integer,Integer,Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetSizes",
    "category": "Method",
    "text": "MatSetSizes\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetType-Tuple{PETSc.PetscMat,Any}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetType",
    "category": "Method",
    "text": "MatSetType\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetValues-Tuple{PETSc.PetscMat,Array{Int64,N},Array{Int64,N},Array{Complex{Float64},N},Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetValues",
    "category": "Method",
    "text": "MatSetValues.  The idx and idy arrays must have PetscInt elements, and the   vals array must have PetscScalar elements.\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatSetValuesBlocked-Tuple{PETSc.PetscMat,Array{Int64,N},Array{Int64,N},Array{Complex{Float64},N},Integer}",
    "page": "Mat Documentation",
    "title": "PETSc.MatSetValuesBlocked",
    "category": "Method",
    "text": "MatSetValues blocked.\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatShellGetContext-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatShellGetContext",
    "category": "Method",
    "text": "MatShellGetContext\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatShellSetOperation-Tuple{PETSc.PetscMat,UInt32,Ptr{Void}}",
    "page": "Mat Documentation",
    "title": "PETSc.MatShellSetOperation",
    "category": "Method",
    "text": "MatShellSetOperation\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatShift-Tuple{PETSc.PetscMat,Number}",
    "page": "Mat Documentation",
    "title": "PETSc.MatShift",
    "category": "Method",
    "text": "MatShift\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatTranspose-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatTranspose",
    "category": "Method",
    "text": "Constructs the transpose of a given matrix.  Can be in place or out of place.\n\nInputs\n\nA: matrix to take the transpose of\n\nKeywords\n\ninplace: whether or not to do the transpose in place\n\nOutputs\n\nA PetscMat object, either the original object A if the transpose was done       in place, or a new matrix object if the transpose was done out of place\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatXAIJSetPreallocation-Tuple{PETSc.PetscMat,Integer,AbstractArray{Int64,1},AbstractArray{Int64,1},AbstractArray{Int64,1},AbstractArray{Int64,1}}",
    "page": "Mat Documentation",
    "title": "PETSc.MatXAIJSetPreallocation",
    "category": "Method",
    "text": "MatXAIJSetPreallocation\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatZeroEntries-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.MatZeroEntries",
    "category": "Method",
    "text": "MatZeroEntries\n\n\n\n"
},

{
    "location": "mat.html#PETSc.PetscDestroy-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "Frees a Petsc object.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "mat.html#PETSc.PetscView",
    "page": "Mat Documentation",
    "title": "PETSc.PetscView",
    "category": "Function",
    "text": "PetscViewer\n\n\n\n"
},

{
    "location": "mat.html#PETSc.SetUp-Tuple{PETSc.PetscMat}",
    "page": "Mat Documentation",
    "title": "PETSc.SetUp",
    "category": "Method",
    "text": "MatSetUp\n\n\n\n"
},

{
    "location": "mat.html#PETSc.MatInfo",
    "page": "Mat Documentation",
    "title": "PETSc.MatInfo",
    "category": "Type",
    "text": "Equivalent to Petsc's MatInfo struct\n\n\n\n"
},

{
    "location": "mat.html#Base.show-Tuple{IO,PETSc.MatInfo}",
    "page": "Mat Documentation",
    "title": "Base.show",
    "category": "Method",
    "text": "show() for a MatInfo\n\n\n\n"
},

{
    "location": "mat.html#Mat-Documentation-1",
    "page": "Mat Documentation",
    "title": "Mat Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page wraps functions in PETSc's Mat API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc]\nPages = [\"mat.jl\"]"
},

{
    "location": "mat_interface.html#",
    "page": "Mat Interface",
    "title": "Mat Interface",
    "category": "page",
    "text": ""
},

{
    "location": "mat_interface.html#PETSc.PetscMat-Tuple{Integer,Integer,Any,MPI.Comm}",
    "page": "Mat Interface",
    "title": "PETSc.PetscMat",
    "category": "Method",
    "text": "Constructor\n\nInputs\n\nmglobal: first dimension global size (or PETSC_DECIDE)\nnglobal: second dimension global size (or PETSC_DECIDE)\nformat: matrix format\ncomm: MPI communicator\n\nKeyword Arguments\n\nmlocal: first dimension local size (or PETSC_DECIDE)\nnlocal: second dimension local size (or PETSC_DECIDE)\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.PetscDestroy-Tuple{AbstractArray{T,2}}",
    "page": "Mat Interface",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "PetscDestroy for AbstractMatrix.  No-op\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.PetscView",
    "page": "Mat Interface",
    "title": "PETSc.PetscView",
    "category": "Function",
    "text": "Print a non-Petsc matrix to a given IO (a Julia IO, not a Petsc IO)\n\nInputs\n\nA: the matrix\nf: an IO, default STDOUT\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.assembly_begin-Tuple{PETSc.PetscMat,Integer}",
    "page": "Mat Interface",
    "title": "PETSc.assembly_begin",
    "category": "Method",
    "text": "Begin matrix assembly for PetscMat.  No-op for Julia matrices\n\nInputs\n\nmat: AbstractMatrix\nflag: type of matrix assembly (see \nMatAssemblyBegin\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.assembly_end-Tuple{PETSc.PetscMat,Integer}",
    "page": "Mat Interface",
    "title": "PETSc.assembly_end",
    "category": "Method",
    "text": "Counterpart of assembly_end\n\nInputs\n\nmat: AbstractMatrix\nflg: type of matrix assembly\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.diagonal_shift!-Tuple{PETSc.PetscMat,Number}",
    "page": "Mat Interface",
    "title": "PETSc.diagonal_shift!",
    "category": "Method",
    "text": "Adds the specified value to the diagonal of the matrix\n\nInputs\n\nA: the matrix\na: the value\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.fill_zero!-Tuple{PETSc.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc.fill_zero!",
    "category": "Method",
    "text": "Fill the matrix with zeros. The sparsity pattern of the matrix (if applicable)   should be defined before this function is called\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.get_values1!-Tuple{AbstractArray{T,2},Array{Int64,N},Array{Int64,N},Array{T,N}}",
    "page": "Mat Interface",
    "title": "PETSc.get_values1!",
    "category": "Method",
    "text": "Method for AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.get_values1!-Tuple{PETSc.PetscMat,Array{Int64,N},Array{Int64,N},Array{Complex{Float64},N}}",
    "page": "Mat Interface",
    "title": "PETSc.get_values1!",
    "category": "Method",
    "text": "Like set_values1!, but retrieves values.  See that function for   the meanings of the arguments. Note that Petsc does   not support getting values for the non-local block of the matrix\n\nInputs\n\nmat: a matrix, can be a Petsc matrix or a julia matrix\n\nInputs/Outputs\n\nidxm\nidxn\nvals\n\nAliasing restrictions: idxm and idxn cannot alias\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.local_indices-Tuple{PETSc.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc.local_indices",
    "category": "Method",
    "text": "Returns the rows owned by this process (1-based)\n\nInputs\n\nA: a matrix\n\nOutputs\n\nrng: a UnitRange\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.set_values1!",
    "page": "Mat Interface",
    "title": "PETSc.set_values1!",
    "category": "Function",
    "text": "1-based indexing for both regular and Pets matrices.   Note that Petsc treats arrays at being row-major, so it is recommened   to set MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE) before using   this function.\n\nInputs\n\nflag: PETSC_INSERT_VALUES or PETSC_ADD_VALUES.  Note that the first one            result in non-deterministic behavior in parallel (in the general            case)\n\nvals: the values, must be length(idxm) x length(idxn)\n\nInputs/Outputs\n\nmat: a matrix, can be a Petsc matrix or a julia matrix\nidxm: the row numbers\nidxn: the column numbers    \n\nNote that idxm and idxn are listed as input/outputs because they may be   modified by this function, however when the function returns they   will have the same values as on entry.  This is needed to accomodate the   fact that Petsc uses 1 based indexing internally.\n\nThis function is optimized for PetscMat and SparseMatrixCSC\n\nAliasing restriction: idxm and idxn cannot alias\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.set_values1!",
    "page": "Mat Interface",
    "title": "PETSc.set_values1!",
    "category": "Function",
    "text": "Method for SparseMatrixCSC\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.set_values1!",
    "page": "Mat Interface",
    "title": "PETSc.set_values1!",
    "category": "Function",
    "text": "Method for AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.size_global-Tuple{PETSc.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc.size_global",
    "category": "Method",
    "text": "Global size of matrix, same as size_local() for serial matrices\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#PETSc.size_local-Tuple{PETSc.PetscMat}",
    "page": "Mat Interface",
    "title": "PETSc.size_local",
    "category": "Method",
    "text": "Size of local part of matrix\n\nInputs\n\nA: AbstractMatrix\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.A_mul_B!-Tuple{PETSc.PetscVec,PETSc.PetscMat,PETSc.PetscVec}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.A_mul_B!",
    "category": "Method",
    "text": "Computes b = A*x\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.At_mul_B!-Tuple{PETSc.PetscVec,PETSc.PetscMat,PETSc.PetscVec}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.At_mul_B!",
    "category": "Method",
    "text": "Computes b = A.'*x\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.norm-Tuple{PETSc.PetscMat,Number}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.norm",
    "category": "Method",
    "text": "Norm for Petsc matrices, 1, 2, and infinity norms supported\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.scale!-Tuple{PETSc.PetscMat,Number}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.scale!",
    "category": "Method",
    "text": "scale! for Petsc matrix\n\n\n\n"
},

{
    "location": "mat_interface.html#Base.LinAlg.vecnorm-Tuple{PETSc.PetscMat}",
    "page": "Mat Interface",
    "title": "Base.LinAlg.vecnorm",
    "category": "Method",
    "text": "Frobenius norm, consistent with Julias interface\n\n\n\n"
},

{
    "location": "mat_interface.html#Mat-Interface-1",
    "page": "Mat Interface",
    "title": "Mat Interface",
    "category": "section",
    "text": "CurrentModule = PETScThis page describes an interface that is efficiently implemented for both Array, SparseMatrixCSC, and Petscmat.  In some cases, Base functions are extended with new methods, in other cases, new functions are defined for the supported matrix types.  Using these functions allows some amount of generic programming.Modules = [PETSc]\nPages = [\"mat_interface.jl\"]"
},

{
    "location": "constants.html#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "constants.html#PETSc.PetscInt",
    "page": "Constants",
    "title": "PETSc.PetscInt",
    "category": "Type",
    "text": "Petsc integer type\n\n\n\n"
},

{
    "location": "constants.html#PETSc.PetscReal",
    "page": "Constants",
    "title": "PETSc.PetscReal",
    "category": "Type",
    "text": "The closest real type to PetscScalar\n\n\n\n"
},

{
    "location": "constants.html#PETSc.PetscScalar",
    "page": "Constants",
    "title": "PETSc.PetscScalar",
    "category": "Type",
    "text": "Element type of Petsc vectors and matrices\n\n\n\n"
},

{
    "location": "constants.html#Constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": "CurrentModule = PETScMany PETSc constants and enums are available in the Julia wrappers.  The full list can be found in src/petsc_constants.jl.  See the PETSc documentation for their meaning. A few of the most important constants are listed here.  Not all of the constants are exported from the PETSc module, to avoid cluttering the users namespace.Modules = [PETSc]\nPages = [\"constants.jl\"]"
},

{
    "location": "ksp.html#",
    "page": "KSP Documentation",
    "title": "KSP Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "ksp.html#PETSc.KSP_NULL",
    "page": "KSP Documentation",
    "title": "PETSc.KSP_NULL",
    "category": "Constant",
    "text": "Null pointer KSP object\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.KSP",
    "page": "KSP Documentation",
    "title": "PETSc.KSP",
    "category": "Type",
    "text": "KSP object\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.GetConvergedReason-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.GetConvergedReason",
    "category": "Method",
    "text": "KSPGetConvergedReason\n\nInputs\n\nksp: KSP object\n\nOutputs\n\nstring containing reason\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.GetInitialGuessNonzero-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.GetInitialGuessNonzero",
    "category": "Method",
    "text": "KSPGetInitialGuessNonzero\n\nInputs\n\nKSP\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.GetResidualNorm-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.GetResidualNorm",
    "category": "Method",
    "text": "KSPGetResidualNorm\n\nInputs\n\nKSP\n\nOutputs\n\nPetscReal\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.GetTolerances-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.GetTolerances",
    "category": "Method",
    "text": "KSPGetTolerances\n\nInputs\n\nKSP\n\nOutputs\n\nrtol\nabstol\ndtol\nmaxits\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.GetType-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.GetType",
    "category": "Method",
    "text": "KSPGetType\n\nInputs\n\nKSP\n\nOutputs\n\nstring containing KSP type\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.KSPSolve-Tuple{PETSc.KSP,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "KSP Documentation",
    "title": "PETSc.KSPSolve",
    "category": "Method",
    "text": "KSPSolve\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.PetscDestroy-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "PetscDestroy for KSP object.  Safe to call multiple times.\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.PetscView",
    "page": "KSP Documentation",
    "title": "PETSc.PetscView",
    "category": "Function",
    "text": "PetscView\n\nInputs\n\nksp: KSP object\nviewer: PetscViewer, defaults to Petsc stdout\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetFromOptions-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.SetFromOptions",
    "category": "Method",
    "text": "SetFromOptions\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetInitialGuessNonzero-Tuple{PETSc.KSP,UInt32}",
    "page": "KSP Documentation",
    "title": "PETSc.SetInitialGuessNonzero",
    "category": "Method",
    "text": "KSPSetInitialGuessNonzero\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetOperators-Tuple{PETSc.KSP,PETSc.PetscMat,PETSc.PetscMat}",
    "page": "KSP Documentation",
    "title": "PETSc.SetOperators",
    "category": "Method",
    "text": "KSPSetOperators\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetTolerances-Tuple{PETSc.KSP,Float64,Float64,Float64,Int64}",
    "page": "KSP Documentation",
    "title": "PETSc.SetTolerances",
    "category": "Method",
    "text": "KSPSetTolerances\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetType-Tuple{PETSc.KSP,ASCIIString}",
    "page": "KSP Documentation",
    "title": "PETSc.SetType",
    "category": "Method",
    "text": "KSPSetType\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetUp-Tuple{PETSc.KSP}",
    "page": "KSP Documentation",
    "title": "PETSc.SetUp",
    "category": "Method",
    "text": "KSPSetUp\n\n\n\n"
},

{
    "location": "ksp.html#PETSc.SetReusePreconditioner-Tuple{PETSc.KSP,UInt32}",
    "page": "KSP Documentation",
    "title": "PETSc.SetReusePreconditioner",
    "category": "Method",
    "text": "KSPSetReusePreconditioner\n\n\n\n"
},

{
    "location": "ksp.html#KSP-Documentation-1",
    "page": "KSP Documentation",
    "title": "KSP Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page wraps functions in PETSc's Vec API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc]\nPages = [\"ksp.jl\"]"
},

{
    "location": "pc.html#",
    "page": "PC Documentation",
    "title": "PC Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "pc.html#PETSc.PC",
    "page": "PC Documentation",
    "title": "PETSc.PC",
    "category": "Type",
    "text": "Petsc PC object\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PC-Tuple{MPI.Comm}",
    "page": "PC Documentation",
    "title": "PETSc.PC",
    "category": "Method",
    "text": "Constructor\n\nInputs\n\ncomm: MPI communicator\n\nOutputs\n\nPC object\n\n\n\n"
},

{
    "location": "pc.html#PETSc.KSPGetPC-Tuple{PETSc.KSP}",
    "page": "PC Documentation",
    "title": "PETSc.KSPGetPC",
    "category": "Method",
    "text": "KSPGetPC\n\nInputs\n\nksp: KSP object\n\nOutput\n\nPC: pc object\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCApply-Tuple{PETSc.PC,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "PC Documentation",
    "title": "PETSc.PCApply",
    "category": "Method",
    "text": "PCApply\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCApplyTranspose-Tuple{PETSc.PC,PETSc.PetscVec,PETSc.PetscVec}",
    "page": "PC Documentation",
    "title": "PETSc.PCApplyTranspose",
    "category": "Method",
    "text": "PCApplyTranspose\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCApplyTransposeExists-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCApplyTransposeExists",
    "category": "Method",
    "text": "PCApplyTransposeExists\n\nInputs\n\npc: PC object\n\nOutputs\n\nBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCBJacobiGetSubKSP-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCBJacobiGetSubKSP",
    "category": "Method",
    "text": "PCBJacobiGetSubKSP\n\nInputs\n\npc: PC object\n\nOutputs\n\nn_local: number of local KSP object\nfirst_local: global number of first KSP object on this block\nksp_arr: array of KSP object\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorGetAllowDiagonalFill-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorGetAllowDiagonalFill",
    "category": "Method",
    "text": "PCFactorGetAllowDiagonalFill\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorGetLevels-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorGetLevels",
    "category": "Method",
    "text": "PCFactorGetLevels\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscInt\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorGetUseInPlace-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorGetUseInPlace",
    "category": "Method",
    "text": "PCFactorGetUseInPlace\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorSetAllowDiagonalFill-Tuple{PETSc.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorSetAllowDiagonalFill",
    "category": "Method",
    "text": "PCFactorSetAllowDiagonalFill\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorSetFill-Tuple{PETSc.PC,Float64}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorSetFill",
    "category": "Method",
    "text": "PCFactorSetFill\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorSetLevels-Tuple{PETSc.PC,Int64}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorSetLevels",
    "category": "Method",
    "text": "PCFactorSetLevels\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCFactorSetUseInPlace-Tuple{PETSc.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc.PCFactorSetUseInPlace",
    "category": "Method",
    "text": "PCFactorSetUseInPlace\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCGetReusePreconditioner-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCGetReusePreconditioner",
    "category": "Method",
    "text": "PCGetReusePreconditioner\n\nInputs\n\npc: PC object\n\nOutputs\n\nPetscBool\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCGetType-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCGetType",
    "category": "Method",
    "text": "PCGetType\n\nInputs\n\npc: PC object\n\nOutputs\n\nstring contianing the PC type\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCJacobiGetType-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCJacobiGetType",
    "category": "Method",
    "text": "PCJacobiGetType\n\nInputs\n\npc: PC object\n\nOutputs\n\nstring containing the type\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCJacobiSetType-Tuple{PETSc.PC,Int32}",
    "page": "PC Documentation",
    "title": "PETSc.PCJacobiSetType",
    "category": "Method",
    "text": "PCJacobiaSetType\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCSetReusePreconditioner-Tuple{PETSc.PC,UInt32}",
    "page": "PC Documentation",
    "title": "PETSc.PCSetReusePreconditioner",
    "category": "Method",
    "text": "PCSetReusePreconditioner\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCSetType-Tuple{PETSc.PC,ASCIIString}",
    "page": "PC Documentation",
    "title": "PETSc.PCSetType",
    "category": "Method",
    "text": "PCSetType\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCSetUp-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCSetUp",
    "category": "Method",
    "text": "PCSetUp\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCShellGetContext-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCShellGetContext",
    "category": "Method",
    "text": "PCShellGetContext\n\nInputs\n\npc: PC object\n\nOutputs\n\nPtr{Void}.  Users shoudl call unsafe_pointer_to_objref() on it to get the                  Julia object back\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCShellSetApply-Tuple{PETSc.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc.PCShellSetApply",
    "category": "Method",
    "text": "PCShellSetApply\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCShellSetApplyTranspose-Tuple{PETSc.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc.PCShellSetApplyTranspose",
    "category": "Method",
    "text": "PCShellSetApplyTranspose\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCShellSetContext-Tuple{PETSc.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc.PCShellSetContext",
    "category": "Method",
    "text": "PCShellSetContext\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCShellSetSetUp-Tuple{PETSc.PC,Ptr{Void}}",
    "page": "PC Documentation",
    "title": "PETSc.PCShellSetSetUp",
    "category": "Method",
    "text": "PCShellSetSetUp\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PetscDestroy-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PetscDestroy",
    "category": "Method",
    "text": "Free a Petsc PC object.  Safe to call multiple times\n\n\n\n"
},

{
    "location": "pc.html#PETSc.KSPSetPC-Tuple{PETSc.KSP,PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.KSPSetPC",
    "category": "Method",
    "text": "KSPSetPC\n\n\n\n"
},

{
    "location": "pc.html#PETSc.PCSetFromOptions-Tuple{PETSc.PC}",
    "page": "PC Documentation",
    "title": "PETSc.PCSetFromOptions",
    "category": "Method",
    "text": "PCSetFromOptions\n\n\n\n"
},

{
    "location": "pc.html#PC-Documentation-1",
    "page": "PC Documentation",
    "title": "PC Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page wraps functions in PETSc's PC API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc]\nPages = [\"pc.jl\"]"
},

{
    "location": "options.html#",
    "page": "Options Documentation",
    "title": "Options Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "options.html#PETSc.PetscClearOptions-Tuple{Dict{K,V}}",
    "page": "Options Documentation",
    "title": "PETSc.PetscClearOptions",
    "category": "Method",
    "text": "Convenience wrapper for using a dictionary to clear options (only the keys   are used).\n\n\n\n"
},

{
    "location": "options.html#PETSc.PetscOptionsClearValue",
    "page": "Options Documentation",
    "title": "PETSc.PetscOptionsClearValue",
    "category": "Function",
    "text": "PetscOptionsClearValue\n\nInputs\n\narg1: the key (string)\narg2: the PetscOptions object, defaults to the global options databse\n\n\n\n"
},

{
    "location": "options.html#PETSc.PetscOptionsSetValue",
    "page": "Options Documentation",
    "title": "PETSc.PetscOptionsSetValue",
    "category": "Function",
    "text": "PetscOptionsSetValue\n\nInputs\n\narg1: the key (string)\narg2: the value (string)\narg3: the PetscOptions object, defaults to the global options database\n\n\n\n"
},

{
    "location": "options.html#PETSc.PetscOptionsView",
    "page": "Options Documentation",
    "title": "PETSc.PetscOptionsView",
    "category": "Function",
    "text": "PetscOptionsView\n\nInputs\n\narg1: a PetscViewer, defaults to Petsc's stdout\narg2: the PetscOptions object, defaults to the global options databse\n\n\n\n"
},

{
    "location": "options.html#PETSc.PetscSetOptions-Tuple{Dict{K,V}}",
    "page": "Options Documentation",
    "title": "PETSc.PetscSetOptions",
    "category": "Method",
    "text": "Convenience wrapper for using a dictionary to set options\n\n\n\n"
},

{
    "location": "options.html#PETSc.PetscOptions",
    "page": "Options Documentation",
    "title": "PETSc.PetscOptions",
    "category": "Type",
    "text": "Typedef of PetscOptions\n\n\n\n"
},

{
    "location": "options.html#Options-Documentation-1",
    "page": "Options Documentation",
    "title": "Options Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page wraps functions in PETSc's options API.  Consult the PETSc documentation for the behavior of the functions.  The documentation on this page only describes differences between the Julia wrappers and the PETSc documentation.Modules = [PETSc]\nPages = [\"options.jl\"]"
},

{
    "location": "error.html#",
    "page": "Error Handling Documentation",
    "title": "Error Handling Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "error.html#PETSc.error_handler-Tuple{MPI.CComm,Int32,Ptr{UInt8},Ptr{UInt8},Int32,Int32,Ptr{UInt8},Ptr{Void}}",
    "page": "Error Handling Documentation",
    "title": "PETSc.error_handler",
    "category": "Method",
    "text": "Error handler registered with PETSc.  This function prints the information   Petsc supplies about the error, a Julia stack trace, and then invokes the   PetscTraceBackErrorHandler to print the Petsc stack trace.\n\n\n\n"
},

{
    "location": "error.html#Error-Handling-Documentation-1",
    "page": "Error Handling Documentation",
    "title": "Error Handling Documentation",
    "category": "section",
    "text": "CurrentModule = PETScThis page describes error handling in the Julia wrappers.Modules = [PETSc]\nPages = [\"error.jl\"]"
},

]}
