
# export names
# do typealiases have to be exported?  I don't think so

export PETSC_NULL, PETSC_IGNORE, PETSC_DECIDE, PETSC_DETERMINE, PETSCDEFAULT, PETSC_COMM_SELF

export PetscInt, PetscScalar, PetscBool, PetscErrorCode, PetscDataType, PetscReal

export PETSC_INSERT_VALUES, PETSC_ADD_VALUES, PETSC_COPY_VALUESA
export PETSC_TRUE, PETSC_FALSE
export NORM_1, NORM_2, NORM_FROBENIUS, NORM_INFINITY, NORM_MAX

export PETSC_MAT_FLUSH_ASSEMBLY, PETSC_MAT_FINAL_ASSEMBLY

export KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGROPPCG, KSPIPECG, KSPCGNE, KSPNASH, KSPSTCG, KSPGLTR, KSPFCG, KSPGMRES, KSPFGMRES, KSPLGMRES, KSPDGMRES, KSPPGMRES, KSPTCQMR, SKPBCGS, KSPIBCGS, KSPFBCGS, KSPFBCGSR, KSPBCGSL, KSPCGS, KSPTFQMR, KSPCR, KSPPIPECR, KSPSLQR, SKPPREONLY, KSPQCQ, KSPBICG, KSPMINRES, KSPSYMMLQ, KSPLCD, KSPPYTHON, KSPCR

export KSP_NORM_DEFAULT, KSP_NORM_NONE, KSP_NORM_PRECONDITIONED, KSP_NORM_UNPRECONGDITIONED, KSP_NORM_NATURAL

export KSP_NORM_MAX



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
  file_path = joinpath(Pkg.dir("PETSc"), "deps/petsc_evars")
  args = open(readdlm, file_path)
  println("args = ", args)
  ENV["PETSC_DIR"] = args[1]
  ENV["PETSC_ARCH"] = args[2]
end

global const PETSC_DIR = ENV["PETSC_DIR"]
global const PETSC_ARCH = ENV["PETSC_ARCH"]


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
typealias Petsc64bitInt Int64
#typealias PetscBLASInt Int32
#typealias PetscScalar Cint
typealias PetscBool Uint32
typealias PetscDataType Cint  # C enums are Int32

typealias PetscLogDouble Cdouble
typealias PetscErrorCode Cint
#=
const PETSC_PRECISION_SINGLE = (UInt32)(4)
const PETSC_PRECISION_DOUBLE = (UInt32)(8)
=#

const PETSC_FALSE = (UInt32)(0)
const PETSC_TRUE = (UInt32)(1)

const PETSC_INT = (Int32)(0)
const PETSC_DOUBLE = (Int32)(1)
const PETSC_COMPLEX = (Int32)(2)
const PETSC_LONG = (Int32)(3)
const PETSC_SHORT = (Int32)(4)
const PETSC_FLOAT = (Int32)(5)
const PETSC_CHAR = (Int32)(6)
const PETSC_BIT_LOGICAL = (Int32)(7)
const PETSC_ENUM = (Int32)(8)
const PETSC_BOOL = (Int32)(9)
const PETSC___FLOAT128 = (Int32)(10)
const PETSC_OBJECT = (Int32)(11)
const PETSC_FUNCTION = (Int32)(12)
const PETSC_STRING = (Int32)(12)



# these functions are used for figuring out the sizes of the datatypes
# use the ones in PETSc.jl for writing programs
function PetscDataTypeFromString_(name::AbstractString)
    ptype = Array(Cint, 1)
    found = Array(PetscBool, 1)
    ccall((:PetscDataTypeFromString,petsc),PetscErrorCode,(Cstring,Ptr{PetscDataType},Ptr{PetscBool}), name, ptype, found)

    return ptype[1], convert(Bool, found[1])
end


function PetscDataTypeGetSize_(dtype::PetscDataType)
    datasize = Array(Csize_t, 1)
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
  println("unknown type of Real")
  println("real_type_enum = ", real_type_enum)
end


if scalar_type_enum == real_type_enum
  scalar_type = "real"
elseif scalar_type_enum == PETSC_COMPLEX
  scalar_type = "complex"
else
  println("unknown type of Scalar")
  println("scalar_type_enum = ", scalar_type_enum)
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
end

println("PetscScalar type = ", scalar_dtype)
println("PetscReal type = ", real_dtype)
println("PetscInt type = ", int_dtype)

typealias PetscScalar scalar_dtype
typealias PetscReal real_dtype
typealias PetscInt int_dtype


# some useful type unions
PetscInt_arr_or_null = Union(AbstractArray{PetscInt}, Ptr{Void})

#=
const PETSC_PI = pi
const PETSC_MAX_INT = 2147483647
const PETSC_MIN_INT = -PETSC_MAX_INT - 1
const PETSC_MAX_REAL = PETSC_MAX_INT  # made up
const PETSC_INFINITY = PETSC_MAX_REAL / 4.0
const PETSC_NINFINITY = -PETSC_INFINITY
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

typealias MPI_Comm MPI.Comm  # use MPI package communicator type
# typealias to size of MPI communicator value
# this is not defined by the MPI C standard, so it might be
# 32 bits (MPICH), or possibly 64 bits
typealias comm_type typeof(MPI.COMM_WORLD.val) 

global const PETSC_INSERT_VALUES = convert(Int32, 1);
global const PETSC_ADD_VALUES    = convert(Int32, 2);
global const PETSC_COPY_VALUES   = convert(Int32, 0);

typealias NormType Int32
global const NORM_1         = convert(Int32, 0);
global const NORM_2         = convert(Int32, 1);
global const NORM_FROBENIUS = convert(Int32 ,2);
global const NORM_INFINITY  = convert(Int32 ,3);
global const NORM_MAX       = NORM_INFINITY;




global const  PETSC_MAT_FLUSH_ASSEMBLY = convert(Int32, 1)
global const  PETSC_MAT_FINAL_ASSEMBLY = convert(Int32, 0)

# vector formats
global const VECSEQ = "seq"
global const VECMPI = "mpi"
global const VECSTANDARD = "standard"
global const VECSHARED = "shared"
global const VECSEQCUSP = "seqcusp"
global const VECMPICUSP = "mpicusp"
global const VECCUSP = "cusp"
global const VECSEQVIENNACL = "seqviennacl"
global const VECMPIVIENNACL = "mpiviennacl"
global const VECVIENNACL = "viennacl"
global const VECNEST = "nest"
global const VECSEQPTHREAD = "seqpthread"
global const VECMPIPTHREAD = "mpipthread"
global const VECPTHREAD = "pthread"

export PetscMatStructure
export DIFFERENT_NONZERO_PATTERN, SUBSET_NONZERO_PATTERN, SAME_NONZERO_PATTERN
# begin enum MatStructure
typealias PetscMatStructure Int32
global const DIFFERENT_NONZERO_PATTERN = (Int32)(0)
global const SUBSET_NONZERO_PATTERN = (Int32)(1)
global const SAME_NONZERO_PATTERN = (Int32)(2)
# end enum MatStructure

export MATSEQAIJ, MATMPIAIJ, MATMPIBAIJ
# matrix types
global const MATSAME = "same"
global const MATMAIJ = "maij"
global const MATSEQMAIJ = "seqmaij"
global const MATMPIMAIJ = "mpimaij"
global const MATIS = "is"
global const MATAIJ = "aij"
global const MATSEQAIJ = "seqaij"
global const MATSEQAIJPTHREAD = "seqaijpthread"
global const MATAIJPTHREAD = "aijpthread"
global const MATMPIAIJ = "mpiaij"
global const MATAIJCRL = "aijcrl"
global const MATSEQAIJCRL = "seqaijcrl"
global const MATMPIAIJCRL = "mpiaijcrl"
global const MATAIJCUSP = "aijcusp"
global const MATSEQAIJCUSP = "seqaijcusp"
global const MATMPIAIJCUSP = "mpiaijcusp"
global const MATAIJCUSPARSE = "aijcusparse"
global const MATSEQAIJCUSPARSE = "seqaijcusparse"
global const MATMPIAIJCUSPARSE = "mpiaijcusparse"
global const MATAIJVIENNACL = "aijviennacl"
global const MATSEQAIJVIENNACL = "seqaijviennacl"
global const MATMPIAIJVIENNACL = "mpiaijviennacl"
global const MATAIJPERM = "aijperm"
global const MATSEQAIJPERM = "seqaijperm"
global const MATMPIAIJPERM = "mpiaijperm"
global const MATSHELL = "shell"
global const MATDENSE = "dense"
global const MATSEQDENSE = "seqdense"
global const MATMPIDENSE = "mpidense"
global const MATELEMENTAL = "elemental"
global const MATBAIJ = "baij"
global const MATSEQBAIJ = "seqbaij"
global const MATMPIBAIJ = "mpibaij"
global const MATMPIADJ = "mpiadj"
global const MATSBAIJ = "sbaij"
global const MATSEQSBAIJ = "seqsbaij"
global const MATMPISBAIJ = "mpisbaij"
global const MATSEQBSTRM = "seqbstrm"
global const MATMPIBSTRM = "mpibstrm"
global const MATBSTRM = "bstrm"
global const MATSEQSBSTRM = "seqsbstrm"
global const MATMPISBSTRM = "mpisbstrm"
global const MATSBSTRM = "sbstrm"
global const MATDAAD = "daad"
global const MATMFFD = "mffd"
global const MATNORMAL = "normal"
global const MATLRC = "lrc"
global const MATSCATTER = "scatter"
global const MATBLOCKMAT = "blockmat"
global const MATCOMPOSITE = "composite"
global const MATFFT = "fft"
global const MATFFTW = "fftw"
global const MATSEQCUFFT = "seqcufft"
global const MATTRANSPOSEMAT = "transpose"
global const MATSCHURCOMPLEMENT = "schurcomplement"
global const MATPYTHON = "python"
global const MATHYPRESTRUCT = "hyprestruct"
global const MATHYPRESSTRUCT = "hypresstruct"
global const MATSUBMATRIX = "submatrix"
global const MATLOCALREF = "localref"
global const MATNEST = "nest"

export MAT_LOCAL, MAT_GLOBAL_MAX, MAT_GLOBAL_SUM
typealias MatInfoType Int32
const MAT_LOCAL = (Int32)(1)
const MAT_GLOBAL_MAX = (Int32)(2)
const MAT_GLOBAL_SUM = (Int32)(3)
#

export MAT_INITIAL_MATRIX, MAT_REUSE_MATRIX, MAT_IGNORE_MATRIX
typealias MatReuse Int32
const MAT_INITIAL_MATRIX = (Int32)(0)
const MAT_REUSE_MATRIX = (Int32)(1)
const MAT_IGNORE_MATRIX = (Int32)(2)
#

export KSPType
# types of KSP solvers
typealias KSPType ASCIIString
global const KSPRICHARDSON = "richardson"
global const KSPCHEBYSHEV = "chebyshev"
global const KSPCG       =  "cg"
global const KSPGROPPCG  =  "groppcg"
global const KSPPIPECG   =  "pipecg"
global const   KSPCGNE   =    "cgne"
global const   KSPNASH   =    "nash"
global const   KSPSTCG   =    "stcg"
global const   KSPGLTR   =    "gltr"
global const KSPFCG      =  "fcg"
global const KSPGMRES    =  "gmres"
global const   KSPFGMRES  =   "fgmres"
global const   KSPLGMRES  =   "lgmres"
global const   KSPDGMRES  =   "dgmres"
global const   KSPPGMRES  =   "pgmres"
global const KSPTCQMR     = "tcqmr"
global const KSPBCGS      = "bcgs"
global const   KSPIBCGS   =   "ibcgs"
global const   KSPFBCGS   =   "fbcgs"
global const   KSPFBCGSR  =   "fbcgsr"
global const   KSPBCGSL   =   "bcgsl"
global const KSPCGS       = "cgs"
global const KSPTFQMR     = "tfqmr"
global const KSPCR        = "cr"
global const KSPPIPECR    = "pipecr"
global const KSPLSQR      = "lsqr"
global const KSPPREONLY   = "preonly"
global const KSPQCG       = "qcg"
global const KSPBICG      = "bicg"
global const KSPMINRES    = "minres"
global const KSPSYMMLQ    = "symmlq"
global const KSPLCD       = "lcd"
global const KSPPYTHON    = "python"
global const KSPGCR       = "gcr"

typealias KSPNormType Cint
global const KSP_NORM_DEFAULT = (Int32)(-1)
global const KSP_NORM_NONE = (Int32)(0)
global const KSP_NORM_PRECONDITIONED = (Int32)(1)
global const KSP_NORM_UNPRECONDITIONED = (Int32)(2)
global const KSP_NORM_NATURAL = (Int32)(3)
# end enum KSPNormType

export KSPConvergedReason

typealias KSPConvergedReason Cint
global const KSP_CONVERGED_RTOL_NORMAL = (Int32)(1)
global const KSP_CONVERGED_ATOL_NORMAL = (Int32)(9)
global const KSP_CONVERGED_RTOL = (Int32)(2)
global const KSP_CONVERGED_ATOL = (Int32)(3)
global const KSP_CONVERGED_ITS = (Int32)(4)
global const KSP_CONVERGED_CG_NEG_CURVE = (Int32)(5)
global const KSP_CONVERGED_CG_CONSTRAINED = (Int32)(6)
global const KSP_CONVERGED_STEP_LENGTH = (Int32)(7)
global const KSP_CONVERGED_HAPPY_BREAKDOWN = (Int32)(8)
global const KSP_DIVERGED_NULL = (Int32)(-2)
global const KSP_DIVERGED_ITS = (Int32)(-3)
global const KSP_DIVERGED_DTOL = (Int32)(-4)
global const KSP_DIVERGED_BREAKDOWN = (Int32)(-5)
global const KSP_DIVERGED_BREAKDOWN_BICG = (Int32)(-6)
global const KSP_DIVERGED_NONSYMMETRIC = (Int32)(-7)
global const KSP_DIVERGED_INDEFINITE_PC = (Int32)(-8)
global const KSP_DIVERGED_NANORINF = (Int32)(-9)
global const KSP_DIVERGED_INDEFINITE_MAT = (Int32)(-10)
global const KSP_DIVERGED_PCSETUP_FAILED = (Int32)(-11)
global const KSP_CONVERGED_ITERATING = (Int32)(0)
#
# map values to string for printing
export KSPConvergedReasonDict
global const KSPConvergedReasonDict = Dict{KSPConvergedReason, ASCIIString} (
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
#


global const KSP_NORM_MAX = KSP_NORM_NATURAL + 1


typealias PCType ASCIIString

global const PCNONE = "none"
global const PCJACOBI = "jacobi"
global const PCSOR = "sor"
global const PCLU = "lu"
global const PCSHELL = "shell"
global const PCBJACOBI = "bjacobi"
global const PCMG = "mg"
global const PCEISENSTAT = "eisenstat"
global const PCILU = "ilu"
global const PCICC = "icc"
global const PCASM = "asm"
global const PCGASM = "gasm"
global const PCKSP = "ksp"
global const PCCOMPOSITE = "composite"
global const PCREDUNDANT = "redundant"
global const PCSPAI = "spai"
global const PCNN = "nn"
global const PCCHOLESKY = "cholesky"
global const PCPBJACOBI = "pbjacobi"
global const PCMAT = "mat"
global const PCHYPRE = "hypre"
global const PCPARMS = "parms"
global const PCFIELDSPLIT = "fieldsplit"
global const PCTFS = "tfs"
global const PCML = "ml"
global const PCGALERKIN = "galerkin"
global const PCEXOTIC = "exotic"
global const PCCP = "cp"
global const PCBFBT = "bfbt"
global const PCLSC = "lsc"
global const PCPYTHON = "python"
global const PCPFMG = "pfmg"
global const PCSYSPFMG = "syspfmg"
global const PCREDISTRIBUTE = "redistribute"
global const PCSVD = "svd"
global const PCGAMG = "gamg"
global const PCSACUSP = "sacusp"
global const PCSACUSPPOLY = "sacusppoly"
global const PCBICGSTABCUSP = "bicgstabcusp"
global const PCAINVCUSP = "ainvcusp"
global const PCBDDC = "bddc"
global const PCKACZMARZ = "kaczmarz"

# begin enum PCSide
typealias PCSide Cint
global const PC_SIDE_DEFAULT = (Int32)(-1)
global const PC_LEFT = (Int32)(0)
global const PC_RIGHT = (Int32)(1)
global const PC_SYMMETRIC = (Int32)(2)
# end enum PCSide

typealias PCJacobiType Int32
global const PC_JACOBI_DIAGONAL = (Int32)(0)
global const PC_JACOBI_ROWMAX = (Int32)(1)
global const PC_JACOBI_ROWSUM = (Int32)(2)
#


typealias PetscErrorType Int32
global const PETSC_ERROR_INITIAL = (Int32)(0)
global const PETSC_ERROR_REPEAT = (Int32)(1)
global const PETSC_ERROR_IN_CXX = (Int32)(2)
#


# begin enum MatOption
typealias MatOption Cint
const MAT_OPTION_MIN = (Int32)(-5)
const MAT_NEW_NONZERO_LOCATION_ERR = (Int32)(-4)
const MAT_UNUSED_NONZERO_LOCATION_ERR = (Int32)(-3)
const MAT_NEW_NONZERO_ALLOCATION_ERR = (Int32)(-2)
const MAT_ROW_ORIENTED = (Int32)(-1)
const MAT_SYMMETRIC = (Int32)(1)
const MAT_STRUCTURALLY_SYMMETRIC = (Int32)(2)
const MAT_NEW_DIAGONALS = (Int32)(3)
const MAT_IGNORE_OFF_PROC_ENTRIES = (Int32)(4)
const MAT_USE_HASH_TABLE = (Int32)(5)
const MAT_KEEP_NONZERO_PATTERN = (Int32)(6)
const MAT_IGNORE_ZERO_ENTRIES = (Int32)(7)
const MAT_USE_INODES = (Int32)(8)
const MAT_HERMITIAN = (Int32)(9)
const MAT_SYMMETRY_ETERNAL = (Int32)(10)
const MAT_DUMMY = (Int32)(11)
const MAT_IGNORE_LOWER_TRIANGULAR = (Int32)(12)
const MAT_ERROR_LOWER_TRIANGULAR = (Int32)(13)
const MAT_GETROW_UPPERTRIANGULAR = (Int32)(14)
const MAT_SPD = (Int32)(15)
const MAT_NO_OFF_PROC_ZERO_ROWS = (Int32)(16)
const MAT_NO_OFF_PROC_ENTRIES = (Int32)(17)
const MAT_NEW_NONZERO_LOCATIONS = (Int32)(18)
const MAT_OPTION_MAX = (Int32)(19)
# end enum MatOption

# begin enum MatOperation
typealias MatOperation Uint32
const MATOP_SET_VALUES = (UInt32)(0)
const MATOP_GET_ROW = (UInt32)(1)
const MATOP_RESTORE_ROW = (UInt32)(2)
const MATOP_MULT = (UInt32)(3)
const MATOP_MULT_ADD = (UInt32)(4)
const MATOP_MULT_TRANSPOSE = (UInt32)(5)
const MATOP_MULT_TRANSPOSE_ADD = (UInt32)(6)
const MATOP_SOLVE = (UInt32)(7)
const MATOP_SOLVE_ADD = (UInt32)(8)
const MATOP_SOLVE_TRANSPOSE = (UInt32)(9)
const MATOP_SOLVE_TRANSPOSE_ADD = (UInt32)(10)
const MATOP_LUFACTOR = (UInt32)(11)
const MATOP_CHOLESKYFACTOR = (UInt32)(12)
const MATOP_SOR = (UInt32)(13)
const MATOP_TRANSPOSE = (UInt32)(14)
const MATOP_GETINFO = (UInt32)(15)
const MATOP_EQUAL = (UInt32)(16)
const MATOP_GET_DIAGONAL = (UInt32)(17)
const MATOP_DIAGONAL_SCALE = (UInt32)(18)
const MATOP_NORM = (UInt32)(19)
const MATOP_ASSEMBLY_BEGIN = (UInt32)(20)
const MATOP_ASSEMBLY_END = (UInt32)(21)
const MATOP_SET_OPTION = (UInt32)(22)
const MATOP_ZERO_ENTRIES = (UInt32)(23)
const MATOP_ZERO_ROWS = (UInt32)(24)
const MATOP_LUFACTOR_SYMBOLIC = (UInt32)(25)
const MATOP_LUFACTOR_NUMERIC = (UInt32)(26)
const MATOP_CHOLESKY_FACTOR_SYMBOLIC = (UInt32)(27)
const MATOP_CHOLESKY_FACTOR_NUMERIC = (UInt32)(28)
const MATOP_SETUP_PREALLOCATION = (UInt32)(29)
const MATOP_ILUFACTOR_SYMBOLIC = (UInt32)(30)
const MATOP_ICCFACTOR_SYMBOLIC = (UInt32)(31)
const MATOP_DUPLICATE = (UInt32)(34)
const MATOP_FORWARD_SOLVE = (UInt32)(35)
const MATOP_BACKWARD_SOLVE = (UInt32)(36)
const MATOP_ILUFACTOR = (UInt32)(37)
const MATOP_ICCFACTOR = (UInt32)(38)
const MATOP_AXPY = (UInt32)(39)
const MATOP_GET_SUBMATRICES = (UInt32)(40)
const MATOP_INCREASE_OVERLAP = (UInt32)(41)
const MATOP_GET_VALUES = (UInt32)(42)
const MATOP_COPY = (UInt32)(43)
const MATOP_GET_ROW_MAX = (UInt32)(44)
const MATOP_SCALE = (UInt32)(45)
const MATOP_SHIFT = (UInt32)(46)
const MATOP_DIAGONAL_SET = (UInt32)(47)
const MATOP_ZERO_ROWS_COLUMNS = (UInt32)(48)
const MATOP_SET_RANDOM = (UInt32)(49)
const MATOP_GET_ROW_IJ = (UInt32)(50)
const MATOP_RESTORE_ROW_IJ = (UInt32)(51)
const MATOP_GET_COLUMN_IJ = (UInt32)(52)
const MATOP_RESTORE_COLUMN_IJ = (UInt32)(53)
const MATOP_FDCOLORING_CREATE = (UInt32)(54)
const MATOP_COLORING_PATCH = (UInt32)(55)
const MATOP_SET_UNFACTORED = (UInt32)(56)
const MATOP_PERMUTE = (UInt32)(57)
const MATOP_SET_VALUES_BLOCKED = (UInt32)(58)
const MATOP_GET_SUBMATRIX = (UInt32)(59)
const MATOP_DESTROY = (UInt32)(60)
const MATOP_VIEW = (UInt32)(61)
const MATOP_CONVERT_FROM = (UInt32)(62)
const MATOP_MATMAT_MULT = (UInt32)(63)
const MATOP_MATMAT_MULT_SYMBOLIC = (UInt32)(64)
const MATOP_MATMAT_MULT_NUMERIC = (UInt32)(65)
const MATOP_SET_LOCAL_TO_GLOBAL_MAP = (UInt32)(66)
const MATOP_SET_VALUES_LOCAL = (UInt32)(67)
const MATOP_ZERO_ROWS_LOCAL = (UInt32)(68)
const MATOP_GET_ROW_MAX_ABS = (UInt32)(69)
const MATOP_GET_ROW_MIN_ABS = (UInt32)(70)
const MATOP_CONVERT = (UInt32)(71)
const MATOP_SET_COLORING = (UInt32)(72)
const MATOP_SET_VALUES_ADIFOR = (UInt32)(74)
const MATOP_FD_COLORING_APPLY = (UInt32)(75)
const MATOP_SET_FROM_OPTIONS = (UInt32)(76)
const MATOP_MULT_CONSTRAINED = (UInt32)(77)
const MATOP_MULT_TRANSPOSE_CONSTRAIN = (UInt32)(78)
const MATOP_FIND_ZERO_DIAGONALS = (UInt32)(79)
const MATOP_MULT_MULTIPLE = (UInt32)(80)
const MATOP_SOLVE_MULTIPLE = (UInt32)(81)
const MATOP_GET_INERTIA = (UInt32)(82)
const MATOP_LOAD = (UInt32)(83)
const MATOP_IS_SYMMETRIC = (UInt32)(84)
const MATOP_IS_HERMITIAN = (UInt32)(85)
const MATOP_IS_STRUCTURALLY_SYMMETRIC = (UInt32)(86)
const MATOP_SET_VALUES_BLOCKEDLOCAL = (UInt32)(87)
const MATOP_GET_VECS = (UInt32)(88)
const MATOP_MAT_MULT = (UInt32)(89)
const MATOP_MAT_MULT_SYMBOLIC = (UInt32)(90)
const MATOP_MAT_MULT_NUMERIC = (UInt32)(91)
const MATOP_PTAP = (UInt32)(92)
const MATOP_PTAP_SYMBOLIC = (UInt32)(93)
const MATOP_PTAP_NUMERIC = (UInt32)(94)
const MATOP_MAT_TRANSPOSE_MULT = (UInt32)(95)
const MATOP_MAT_TRANSPOSE_MULT_SYMBO = (UInt32)(96)
const MATOP_MAT_TRANSPOSE_MULT_NUMER = (UInt32)(97)
const MATOP_CONJUGATE = (UInt32)(102)
const MATOP_SET_VALUES_ROW = (UInt32)(104)
const MATOP_REAL_PART = (UInt32)(105)
const MATOP_IMAGINARY_PART = (UInt32)(106)
const MATOP_GET_ROW_UPPER_TRIANGULAR = (UInt32)(107)
const MATOP_RESTORE_ROW_UPPER_TRIANG = (UInt32)(108)
const MATOP_MAT_SOLVE = (UInt32)(109)
const MATOP_GET_REDUNDANT_MATRIX = (UInt32)(110)
const MATOP_GET_ROW_MIN = (UInt32)(111)
const MATOP_GET_COLUMN_VECTOR = (UInt32)(112)
const MATOP_MISSING_DIAGONAL = (UInt32)(113)
const MATOP_GET_SEQ_NONZERO_STRUCTUR = (UInt32)(114)
const MATOP_CREATE = (UInt32)(115)
const MATOP_GET_GHOSTS = (UInt32)(116)
const MATOP_GET_LOCAL_SUB_MATRIX = (UInt32)(117)
const MATOP_RESTORE_LOCALSUB_MATRIX = (UInt32)(118)
const MATOP_MULT_DIAGONAL_BLOCK = (UInt32)(119)
const MATOP_HERMITIAN_TRANSPOSE = (UInt32)(120)
const MATOP_MULT_HERMITIAN_TRANSPOSE = (UInt32)(121)
const MATOP_MULT_HERMITIAN_TRANS_ADD = (UInt32)(122)
const MATOP_GET_MULTI_PROC_BLOCK = (UInt32)(123)
const MATOP_FIND_NONZERO_ROWS = (UInt32)(124)
const MATOP_GET_COLUMN_NORMS = (UInt32)(125)
const MATOP_INVERT_BLOCK_DIAGONAL = (UInt32)(126)
const MATOP_GET_SUB_MATRICES_PARALLE = (UInt32)(128)
const MATOP_SET_VALUES_BATCH = (UInt32)(129)
const MATOP_TRANSPOSE_MAT_MULT = (UInt32)(130)
const MATOP_TRANSPOSE_MAT_MULT_SYMBO = (UInt32)(131)
const MATOP_TRANSPOSE_MAT_MULT_NUMER = (UInt32)(132)
const MATOP_TRANSPOSE_COLORING_CREAT = (UInt32)(133)
const MATOP_TRANS_COLORING_APPLY_SPT = (UInt32)(134)
const MATOP_TRANS_COLORING_APPLY_DEN = (UInt32)(135)
const MATOP_RART = (UInt32)(136)
const MATOP_RART_SYMBOLIC = (UInt32)(137)
const MATOP_RART_NUMERIC = (UInt32)(138)
const MATOP_SET_BLOCK_SIZES = (UInt32)(139)
const MATOP_AYPX = (UInt32)(140)
const MATOP_RESIDUAL = (UInt32)(141)
const MATOP_FDCOLORING_SETUP = (UInt32)(142)
const MATOP_MPICONCATENATESEQ = (UInt32)(144)
# end enum MatOperation


