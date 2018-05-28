
# export names
# do typealiases have to be exported?  I don't think so

global const comm_type = MPI.CComm

export PETSC_NULL, PETSC_IGNORE, PETSC_DECIDE, PETSC_DETERMINE, PETSCDEFAULT, PETSC_COMM_SELF

export PetscInt, PetscScalar, PetscBool, PetscErrorCode, PetscDataType, PetscReal

export INSERT_VALUES, ADD_VALUES, PETSC_COPY_VALUESA
export PETSC_TRUE, PETSC_FALSE
export NORM_1, NORM_2, NORM_FROBENIUS, NORM_INFINITY, NORM_MAX

export MAT_FLUSH_ASSEMBLY, MAT_FINAL_ASSEMBLY

export KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGROPPCG, KSPIPECG, KSPCGNE, KSPNASH, KSPSTCG, KSPGLTR, KSPFCG, KSPGMRES, KSPFGMRES, KSPLGMRES, KSPDGMRES, KSPPGMRES, KSPTCQMR, SKPBCGS, KSPIBCGS, KSPFBCGS, KSPFBCGSR, KSPBCGSL, KSPCGS, KSPTFQMR, KSPCR, KSPPIPECR, KSPSLQR, SKPPREONLY, KSPQCQ, KSPBICG, KSPMINRES, KSPSYMMLQ, KSPLCD, KSPPYTHON, KSPCR

export KSP_NORM_DEFAULT, KSP_NORM_NONE, KSP_NORM_PRECONDITIONED, KSP_NORM_UNPRECONGDITIONED, KSP_NORM_NATURAL

export KSP_NORM_MAX

export PetscIntNullArray, PetscScalarNullArray


export VECSEQ
export VECMPI
export VECSTANDARD
export VECSHARED
export VECSEQCUSP
export VECMPICUSP
export VECCUSP
export VECSEQVIENNACL
export VECMPIVIENNACL
export VECVIENNACL
export VECNEST
export VECSEQPTHREAD
export VECMPIPTHREAD
export VECPTHREAD



# decide which version of Petsc to use
# if no PETSC_DIR or PETSC_ARCH defined, use the one (hopefully) built 
# during package installation
if !haskey(ENV, "PETSC_DIR") && !haskey(ENV, "PETSC_ARCH")
  file_path = joinpath(Pkg.dir("PETSc2"), "deps/petsc_evars")
  args = open(readdlm, file_path)
  ENV["PETSC_DIR"] = args[1]
  ENV["PETSC_ARCH"] = args[2]
end

global const PETSC_DIR = ENV["PETSC_DIR"]
global const PETSC_ARCH = ENV["PETSC_ARCH"]

function finalizeMPI()
  if MPI.Initialized()
    MPI.Finalize()
  end
end
if !MPI.Initialized()
  MPI.Init()
  atexit(finalizeMPI)
end
  
myrank = MPI.Comm_rank(MPI.COMM_WORLD)

#PETSC_DIR = readall(`echo $PETSC_DIR`)
#PETSC_ARCH = readall(`echo $PETSC_ARCH`)

#PETSC_DIR = "/home/jared/build/petsc-3.6.0"
#PETSC_ARCH = "arch-linux2-c-debug"

#PETSC_DIR = getenv("PETSC_DIR");
#PETSC_ARCH = getenv("PETSC_ARCH");
#=
if (length(PETSC_DIR) == 0)
  disp("Must have environmental variable PETSC_DIR set")
end
if (length(PETSC_ARCH) == 0)
  disp("Must have environmental variable PETSC_ARCH set")
end
=#

global const libpetsclocation = string(PETSC_DIR, "/", PETSC_ARCH, "/lib/", "libpetsc")
global const petsc = libpetsclocation # for compatability with auto generated wrappers
global const libpetsc = Libdl.dlopen(libpetsclocation)



# definitions of Petsc types
# a way to automatically generate these would be preferable
# define all the data types that are known a priori
global const Petsc64bitInt = Int64
#global const PetscBLASInt = Int32
#global const PetscScalar = Cint
global const PetscBool = UInt32
global const PetscDataType = Cint  # C enums are Int32

global const PetscLogDouble = Cdouble
global const PetscErrorCode = Cint
#=
global const PETSC_PRECISION_SINGLE = (UInt32)(4)
global const PETSC_PRECISION_DOUBLE = (UInt32)(8)
=#

global const PetscViewer = Ptr{Void}

global const PETSC_FALSE = (UInt32)(0)
global const PETSC_TRUE = (UInt32)(1)

global const PETSC_INT = (Int32)(0)
global const PETSC_DOUBLE = (Int32)(1)
global const PETSC_COMPLEX = (Int32)(2)
global const PETSC_LONG = (Int32)(3)
global const PETSC_SHORT = (Int32)(4)
global const PETSC_FLOAT = (Int32)(5)
global const PETSC_CHAR = (Int32)(6)
global const PETSC_BIT_LOGICAL = (Int32)(7)
global const PETSC_ENUM = (Int32)(8)
global const PETSC_BOOL = (Int32)(9)
global const PETSC___FLOAT128 = (Int32)(10)
global const PETSC_OBJECT = (Int32)(11)
global const PETSC_FUNCTION = (Int32)(12)
global const PETSC_STRING = (Int32)(12)



# these functions are used for figuring out the sizes of the datatypes
# use the ones in PETSc.jl for writing programs
function PetscDataTypeFromString_(name::AbstractString)
    ptype = Array{Cint}(1)
    found = Array{PetscBool}(1)
    ccall((:PetscDataTypeFromString,petsc),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}), name, ptype, found)

    return ptype[1], convert(Bool, found[1])
end


function PetscDataTypeGetSize_(dtype::PetscDataType)
    datasize = Array{Csize_t}(1)
    ccall((:PetscDataTypeGetSize,petsc),PetscErrorCode,(PetscDataType,Ptr{Csize_t}), dtype, datasize)

    return datasize[1]
end


# define types that depend on the options Petsc was compiled with
(real_type_enum, found_real) = PetscDataTypeFromString_("Real")
(scalar_type_enum, found_scalar) = PetscDataTypeFromString_("Scalar")
int_size = PetscDataTypeGetSize_(PETSC_INT)

# confirm a value was found
@assert(found_real)
@assert(found_scalar)

# figure out what the values mean
if real_type_enum == PETSC_DOUBLE
  precision = "double"
elseif real_type_enum == PETSC_FLOAT
  precision = "single"
else
  if myrank == 0
    println("unknown type of Real")
    println("real_type_enum = ", real_type_enum)
  end
  throw(ErrorException("unsupported PetscReal ty pe"))
end


if scalar_type_enum == real_type_enum
  scalar_type = "real"
elseif scalar_type_enum == PETSC_COMPLEX
  scalar_type = "complex"
else
  if myrank == 0
    println("unknown type of Scalar")
    println("scalar_type_enum = ", scalar_type_enum)
  end
  throw(ErrorException("unsupported PetscScalar type"))
end

# figure out types

if scalar_type == "real"
  # single or double precision real
  if precision == "single"
    scalar_dtype = Float32
    real_dtype = Float32
  else
    scalar_dtype = Float64
    real_dtype = Float64
  end
else  # scalar_type = complex
  if precision == "single"
    scalar_dtype = Complex64
    real_dtype = Float32
  else
    scalar_dtype = Complex128
    real_dtype = Float64
  end
end

if int_size == 4
   int_dtype = Int32
elseif int_size == 8
   int_dtype = Int64
else
  println("unknown Int size")
  println("int_size = ", int_size)
  throw(ErrorException("unsupported integer size"))
end


if myrank == 0
  println("PetscScalar type = ", scalar_dtype)
  println("PetscReal type = ", real_dtype)
  println("PetscInt type = ", int_dtype)
end

"""
  Element type of Petsc vectors and matrices
"""
global const PetscScalar = scalar_dtype

"""
  The closest real type to `PetscScalar`
"""
global const PetscReal = real_dtype

"""
  Petsc integer type
"""
global const PetscInt = int_dtype


# some useful type unions
PetscInt_arr_or_null = Union{AbstractArray{PetscInt}, Ptr{Void}}

global const PetscIntNullArray = unsafe_wrap(Vector{PetscInt}, Ptr{PetscInt}(C_NULL), 0)
global const PetscScalarNullArray = unsafe_wrap(Vector{PetscScalar}, Ptr{PetscScalar}(C_NULL), 0)

#=
global const PETSC_PI = pi
global const PETSC_MAX_INT = 2147483647
global const PETSC_MIN_INT = -PETSC_MAX_INT - 1
global const PETSC_MAX_REAL = PETSC_MAX_INT  # made up
global const PETSC_INFINITY = PETSC_MAX_REAL / 4.0
global const PETSC_NINFINITY = -PETSC_INFINITY
#const PassiveReal = Cint  # what is this?
=#

# some useful constants
#const PassiveScalar = PetscScalar
#const MPIU_MATSCALAR = MPIU_SCALAR
#const MPIU_2INT = Cint
global const PETSC_NULL = C_NULL
global const PETSC_IGNORE = C_NULL
global const PETSC_DECIDE = convert(Int32, -1)
global const PETSC_DETERMINE = PETSC_DECIDE
global const PETSC_DEFAULT = convert(Int32, -2)
global const PETSC_COMM_SELF = MPI.COMM_SELF

global const MPI_Comm = MPI.Comm  # use MPI package communicator type
# typealias to size of MPI communicator value
# this is not defined by the MPI C standard, so it might be
# 32 bits (MPICH), or possibly 64 bits
global const comm_type = MPI.CComm

export PetscMatStructure
export DIFFERENT_NONZERO_PATTERN, SUBSET_NONZERO_PATTERN, SAME_NONZERO_PATTERN
export MATSEQAIJ, MATMPIAIJ, MATMPIBAIJ
export MAT_LOCAL, MAT_GLOBAL_MAX, MAT_GLOBAL_SUM
export MAT_INITIAL_MATRIX, MAT_REUSE_MATRIX, MAT_IGNORE_MATRIX
export KSPType
export KSPConvergedReason

# get the auto generated definitions
include("petsc_constants_gen.jl")

# map values to string for printing
export KSPConvergedReasonDict
global const KSPConvergedReasonDict = Dict{KSPConvergedReason, String}(
 KSP_CONVERGED_RTOL_NORMAL => "Converged: RTol normal",
 KSP_CONVERGED_ATOL_NORMAL => "Converged: ATol normal",
 KSP_CONVERGED_RTOL => "Converged: RTol",
 KSP_CONVERGED_ATOL => "converged: ATol",
 KSP_CONVERGED_ITS => "Converged: ITS",
 KSP_CONVERGED_CG_NEG_CURVE => "Converged: CG Negative Curvature",
 KSP_CONVERGED_CG_CONSTRAINED => "Converged: CG Constrained",
 KSP_CONVERGED_STEP_LENGTH => "Converged: step length",
 KSP_CONVERGED_HAPPY_BREAKDOWN => "Converged: happy breakdown",
 KSP_DIVERGED_NULL => "Diverged: Null",
 KSP_DIVERGED_ITS => "Diverged: iteration max",
 KSP_DIVERGED_DTOL => "Diverged divergence tolerance",
 KSP_DIVERGED_BREAKDOWN => "Diverged: Krylov breakdown",
 KSP_DIVERGED_BREAKDOWN_BICG => "Diverged: BICG breakdown",
 KSP_DIVERGED_NONSYMMETRIC => "Diverged: non symmetric",
 KSP_DIVERGED_INDEFINITE_PC => "Diverged: indefinitate PC",
 KSP_DIVERGED_NANORINF => "Diverged nan or inf",
 KSP_DIVERGED_INDEFINITE_MAT => "Diverged: indefinate matrix",
 KSP_DIVERGED_PCSETUP_FAILED => "Diverged: PC setup failed",
 KSP_CONVERGED_ITERATING => "Still running: please be patient",
 )

