typealias PetscErrorType Int32
global const PETSC_ERROR_INITIAL = Int32( 0 )
global const PETSC_ERROR_REPEAT = Int32( 1 )
global const PETSC_ERROR_IN_CXX = Int32( 2 )

typealias KSPConvergedReason Cint
global const KSP_CONVERGED_RTOL_NORMAL = Cint( 1 )
global const KSP_CONVERGED_ATOL_NORMAL = Cint( 9 )
global const KSP_CONVERGED_RTOL = Cint( 2 )
global const KSP_CONVERGED_ATOL = Cint( 3 )
global const KSP_CONVERGED_ITS = Cint( 4 )
global const KSP_CONVERGED_CG_NEG_CURVE = Cint( 5 )
global const KSP_CONVERGED_CG_CONSTRAINED = Cint( 6 )
global const KSP_CONVERGED_STEP_LENGTH = Cint( 7 )
global const KSP_CONVERGED_HAPPY_BREAKDOWN = Cint( 8 )
global const KSP_DIVERGED_NULL = Cint( -2 )
global const KSP_DIVERGED_ITS = Cint( -3 )
global const KSP_DIVERGED_DTOL = Cint( -4 )
global const KSP_DIVERGED_BREAKDOWN = Cint( -5 )
global const KSP_DIVERGED_BREAKDOWN_BICG = Cint( -6 )
global const KSP_DIVERGED_NONSYMMETRIC = Cint( -7 )
global const KSP_DIVERGED_INDEFINITE_PC = Cint( -8 )
global const KSP_DIVERGED_NANORINF = Cint( -9 )
global const KSP_DIVERGED_INDEFINITE_MAT = Cint( -10 )
global const KSP_DIVERGED_PCSETUP_FAILED = Cint( -11 )
global const KSP_CONVERGED_ITERATING = Cint( 0 )

typealias MatOption Cint
global const MAT_OPTION_MIN = Cint( -3 )
global const MAT_NEW_NONZERO_LOCATION_ERR = Cint( 11 )
global const MAT_UNUSED_NONZERO_LOCATION_ERR = Cint( -2 )
global const MAT_NEW_NONZERO_ALLOCATION_ERR = Cint( 19 )
global const MAT_ROW_ORIENTED = Cint( -1 )
global const MAT_SYMMETRIC = Cint( 1 )
global const MAT_STRUCTURALLY_SYMMETRIC = Cint( 2 )
global const MAT_NEW_DIAGONALS = Cint( 3 )
global const MAT_IGNORE_OFF_PROC_ENTRIES = Cint( 4 )
global const MAT_USE_HASH_TABLE = Cint( 5 )
global const MAT_KEEP_NONZERO_PATTERN = Cint( 6 )
global const MAT_IGNORE_ZERO_ENTRIES = Cint( 7 )
global const MAT_USE_INODES = Cint( 8 )
global const MAT_HERMITIAN = Cint( 9 )
global const MAT_SYMMETRY_ETERNAL = Cint( 10 )
global const MAT_IGNORE_LOWER_TRIANGULAR = Cint( 12 )
global const MAT_ERROR_LOWER_TRIANGULAR = Cint( 13 )
global const MAT_GETROW_UPPERTRIANGULAR = Cint( 14 )
global const MAT_SPD = Cint( 15 )
global const MAT_NO_OFF_PROC_ZERO_ROWS = Cint( 16 )
global const MAT_NO_OFF_PROC_ENTRIES = Cint( 17 )
global const MAT_NEW_NONZERO_LOCATIONS = Cint( 18 )
global const MAT_OPTION_MAX = Cint( 21 )

typealias MatInfoType Int32
global const MAT_LOCAL = Int32( 1 )
global const MAT_GLOBAL_MAX = Int32( 2 )
global const MAT_GLOBAL_SUM = Int32( 3 )

typealias KSPNormType Cint
global const KSP_NORM_DEFAULT = Cint( -1 )
global const KSP_NORM_NONE = Cint( 0 )
global const KSP_NORM_PRECONDITIONED = Cint( 1 )
global const KSP_NORM_UNPRECONDITIONED = Cint( 2 )
global const KSP_NORM_NATURAL = Cint( 3 )

typealias MatOperation UInt32
global const MATOP_SET_VALUES = UInt32( 0 )
global const MATOP_GET_ROW = UInt32( 1 )
global const MATOP_RESTORE_ROW = UInt32( 2 )
global const MATOP_MULT = UInt32( 3 )
global const MATOP_MULT_ADD = UInt32( 4 )
global const MATOP_MULT_TRANSPOSE = UInt32( 5 )
global const MATOP_MULT_TRANSPOSE_ADD = UInt32( 6 )
global const MATOP_SOLVE = UInt32( 7 )
global const MATOP_SOLVE_ADD = UInt32( 8 )
global const MATOP_SOLVE_TRANSPOSE = UInt32( 9 )
global const MATOP_SOLVE_TRANSPOSE_ADD = UInt32( 10 )
global const MATOP_LUFACTOR = UInt32( 11 )
global const MATOP_CHOLESKYFACTOR = UInt32( 12 )
global const MATOP_SOR = UInt32( 13 )
global const MATOP_TRANSPOSE = UInt32( 14 )
global const MATOP_GETINFO = UInt32( 15 )
global const MATOP_EQUAL = UInt32( 16 )
global const MATOP_GET_DIAGONAL = UInt32( 17 )
global const MATOP_DIAGONAL_SCALE = UInt32( 18 )
global const MATOP_NORM = UInt32( 19 )
global const MATOP_ASSEMBLY_BEGIN = UInt32( 20 )
global const MATOP_ASSEMBLY_END = UInt32( 21 )
global const MATOP_SET_OPTION = UInt32( 22 )
global const MATOP_ZERO_ENTRIES = UInt32( 23 )
global const MATOP_ZERO_ROWS = UInt32( 24 )
global const MATOP_LUFACTOR_SYMBOLIC = UInt32( 25 )
global const MATOP_LUFACTOR_NUMERIC = UInt32( 26 )
global const MATOP_CHOLESKY_FACTOR_SYMBOLIC = UInt32( 27 )
global const MATOP_CHOLESKY_FACTOR_NUMERIC = UInt32( 28 )
global const MATOP_SETUP_PREALLOCATION = UInt32( 29 )
global const MATOP_ILUFACTOR_SYMBOLIC = UInt32( 30 )
global const MATOP_ICCFACTOR_SYMBOLIC = UInt32( 31 )
global const MATOP_DUPLICATE = UInt32( 34 )
global const MATOP_FORWARD_SOLVE = UInt32( 35 )
global const MATOP_BACKWARD_SOLVE = UInt32( 36 )
global const MATOP_ILUFACTOR = UInt32( 37 )
global const MATOP_ICCFACTOR = UInt32( 38 )
global const MATOP_AXPY = UInt32( 39 )
global const MATOP_GET_SUBMATRICES = UInt32( 40 )
global const MATOP_INCREASE_OVERLAP = UInt32( 41 )
global const MATOP_GET_VALUES = UInt32( 42 )
global const MATOP_COPY = UInt32( 43 )
global const MATOP_GET_ROW_MAX = UInt32( 44 )
global const MATOP_SCALE = UInt32( 45 )
global const MATOP_SHIFT = UInt32( 46 )
global const MATOP_DIAGONAL_SET = UInt32( 47 )
global const MATOP_ZERO_ROWS_COLUMNS = UInt32( 48 )
global const MATOP_SET_RANDOM = UInt32( 49 )
global const MATOP_GET_ROW_IJ = UInt32( 50 )
global const MATOP_RESTORE_ROW_IJ = UInt32( 51 )
global const MATOP_GET_COLUMN_IJ = UInt32( 52 )
global const MATOP_RESTORE_COLUMN_IJ = UInt32( 53 )
global const MATOP_FDCOLORING_CREATE = UInt32( 54 )
global const MATOP_COLORING_PATCH = UInt32( 55 )
global const MATOP_SET_UNFACTORED = UInt32( 56 )
global const MATOP_PERMUTE = UInt32( 57 )
global const MATOP_SET_VALUES_BLOCKED = UInt32( 58 )
global const MATOP_GET_SUBMATRIX = UInt32( 59 )
global const MATOP_DESTROY = UInt32( 60 )
global const MATOP_VIEW = UInt32( 61 )
global const MATOP_CONVERT_FROM = UInt32( 62 )
global const MATOP_MATMAT_MULT = UInt32( 63 )
global const MATOP_MATMAT_MULT_SYMBOLIC = UInt32( 64 )
global const MATOP_MATMAT_MULT_NUMERIC = UInt32( 65 )
global const MATOP_SET_LOCAL_TO_GLOBAL_MAP = UInt32( 66 )
global const MATOP_SET_VALUES_LOCAL = UInt32( 67 )
global const MATOP_ZERO_ROWS_LOCAL = UInt32( 68 )
global const MATOP_GET_ROW_MAX_ABS = UInt32( 69 )
global const MATOP_GET_ROW_MIN_ABS = UInt32( 70 )
global const MATOP_CONVERT = UInt32( 71 )
global const MATOP_SET_COLORING = UInt32( 72 )
global const MATOP_SET_VALUES_ADIFOR = UInt32( 74 )
global const MATOP_FD_COLORING_APPLY = UInt32( 75 )
global const MATOP_SET_FROM_OPTIONS = UInt32( 76 )
global const MATOP_MULT_CONSTRAINED = UInt32( 77 )
global const MATOP_MULT_TRANSPOSE_CONSTRAIN = UInt32( 78 )
global const MATOP_FIND_ZERO_DIAGONALS = UInt32( 79 )
global const MATOP_MULT_MULTIPLE = UInt32( 80 )
global const MATOP_SOLVE_MULTIPLE = UInt32( 81 )
global const MATOP_GET_INERTIA = UInt32( 82 )
global const MATOP_LOAD = UInt32( 83 )
global const MATOP_IS_SYMMETRIC = UInt32( 84 )
global const MATOP_IS_HERMITIAN = UInt32( 85 )
global const MATOP_IS_STRUCTURALLY_SYMMETRIC = UInt32( 86 )
global const MATOP_SET_VALUES_BLOCKEDLOCAL = UInt32( 87 )
global const MATOP_GET_VECS = UInt32( 88 )
global const MATOP_MAT_MULT = UInt32( 89 )
global const MATOP_MAT_MULT_SYMBOLIC = UInt32( 90 )
global const MATOP_MAT_MULT_NUMERIC = UInt32( 91 )
global const MATOP_PTAP = UInt32( 92 )
global const MATOP_PTAP_SYMBOLIC = UInt32( 93 )
global const MATOP_PTAP_NUMERIC = UInt32( 94 )
global const MATOP_MAT_TRANSPOSE_MULT = UInt32( 95 )
global const MATOP_MAT_TRANSPOSE_MULT_SYMBO = UInt32( 96 )
global const MATOP_MAT_TRANSPOSE_MULT_NUMER = UInt32( 97 )
global const MATOP_CONJUGATE = UInt32( 102 )
global const MATOP_SET_VALUES_ROW = UInt32( 104 )
global const MATOP_REAL_PART = UInt32( 105 )
global const MATOP_IMAGINARY_PART = UInt32( 106 )
global const MATOP_GET_ROW_UPPER_TRIANGULAR = UInt32( 107 )
global const MATOP_RESTORE_ROW_UPPER_TRIANG = UInt32( 108 )
global const MATOP_MAT_SOLVE = UInt32( 109 )
global const MATOP_GET_REDUNDANT_MATRIX = UInt32( 110 )
global const MATOP_GET_ROW_MIN = UInt32( 111 )
global const MATOP_GET_COLUMN_VECTOR = UInt32( 112 )
global const MATOP_MISSING_DIAGONAL = UInt32( 113 )
global const MATOP_GET_SEQ_NONZERO_STRUCTUR = UInt32( 114 )
global const MATOP_CREATE = UInt32( 115 )
global const MATOP_GET_GHOSTS = UInt32( 116 )
global const MATOP_GET_LOCAL_SUB_MATRIX = UInt32( 117 )
global const MATOP_RESTORE_LOCALSUB_MATRIX = UInt32( 118 )
global const MATOP_MULT_DIAGONAL_BLOCK = UInt32( 119 )
global const MATOP_HERMITIAN_TRANSPOSE = UInt32( 120 )
global const MATOP_MULT_HERMITIAN_TRANSPOSE = UInt32( 121 )
global const MATOP_MULT_HERMITIAN_TRANS_ADD = UInt32( 122 )
global const MATOP_GET_MULTI_PROC_BLOCK = UInt32( 123 )
global const MATOP_FIND_NONZERO_ROWS = UInt32( 124 )
global const MATOP_GET_COLUMN_NORMS = UInt32( 125 )
global const MATOP_INVERT_BLOCK_DIAGONAL = UInt32( 126 )
global const MATOP_GET_SUB_MATRICES_PARALLE = UInt32( 128 )
global const MATOP_SET_VALUES_BATCH = UInt32( 129 )
global const MATOP_TRANSPOSE_MAT_MULT = UInt32( 130 )
global const MATOP_TRANSPOSE_MAT_MULT_SYMBO = UInt32( 131 )
global const MATOP_TRANSPOSE_MAT_MULT_NUMER = UInt32( 132 )
global const MATOP_TRANSPOSE_COLORING_CREAT = UInt32( 133 )
global const MATOP_TRANS_COLORING_APPLY_SPT = UInt32( 134 )
global const MATOP_TRANS_COLORING_APPLY_DEN = UInt32( 135 )
global const MATOP_RART = UInt32( 136 )
global const MATOP_RART_SYMBOLIC = UInt32( 137 )
global const MATOP_RART_NUMERIC = UInt32( 138 )
global const MATOP_SET_BLOCK_SIZES = UInt32( 139 )
global const MATOP_AYPX = UInt32( 140 )
global const MATOP_RESIDUAL = UInt32( 141 )
global const MATOP_FDCOLORING_SETUP = UInt32( 142 )
global const MATOP_MPICONCATENATESEQ = UInt32( 144 )

typealias InsertMode Cint
global const INSERT_VALUES = Cint( 1 )
global const ADD_VALUES = Cint( 2 )

typealias KSPType String
global const KSPRICHARDSON = String( "richardson" )
global const KSPCHEBYSHEV = String( "chebyshev" )
global const KSPCG = String( "cg" )
global const KSPGROPPCG = String( "groppcg" )
global const KSPPIPECG = String( "pipecg" )
global const KSPCGNE = String( "cgne" )
global const KSPNASH = String( "nash" )
global const KSPSTCG = String( "stcg" )
global const KSPGLTR = String( "gltr" )
global const KSPFCG = String( "fcg" )
global const KSPGMRES = String( "gmres" )
global const KSPFGMRES = String( "fgmres" )
global const KSPLGMRES = String( "lgmres" )
global const KSPDGMRES = String( "dgmres" )
global const KSPPGMRES = String( "pgmres" )
global const KSPTCQMR = String( "tcqmr" )
global const KSPBCGS = String( "bcgs" )
global const KSPIBCGS = String( "ibcgs" )
global const KSPFBCGS = String( "fbcgs" )
global const KSPFBCGSR = String( "fbcgsr" )
global const KSPBCGSL = String( "bcgsl" )
global const KSPCGS = String( "cgs" )
global const KSPTFQMR = String( "tfqmr" )
global const KSPCR = String( "cr" )
global const KSPPIPECR = String( "pipecr" )
global const KSPLSQR = String( "lsqr" )
global const KSPPREONLY = String( "preonly" )
global const KSPQCG = String( "qcg" )
global const KSPBICG = String( "bicg" )
global const KSPMINRES = String( "minres" )
global const KSPSYMMLQ = String( "symmlq" )
global const KSPLCD = String( "lcd" )
global const KSPPYTHON = String( "python" )
global const KSPGCR = String( "gcr" )

typealias PCJacobiType Int32
global const PC_JACOBI_DIAGONAL = Int32( 0 )
global const PC_JACOBI_ROWMAX = Int32( 1 )
global const PC_JACOBI_ROWSUM = Int32( 2 )

typealias PCType String
global const PCNONE = String( "none" )
global const PCJACOBI = String( "jacobi" )
global const PCSOR = String( "sor" )
global const PCLU = String( "lu" )
global const PCSHELL = String( "shell" )
global const PCBJACOBI = String( "bjacobi" )
global const PCMG = String( "mg" )
global const PCEISENSTAT = String( "eisenstat" )
global const PCILU = String( "ilu" )
global const PCICC = String( "icc" )
global const PCASM = String( "asm" )
global const PCGASM = String( "gasm" )
global const PCKSP = String( "ksp" )
global const PCCOMPOSITE = String( "composite" )
global const PCREDUNDANT = String( "redundant" )
global const PCSPAI = String( "spai" )
global const PCNN = String( "nn" )
global const PCCHOLESKY = String( "cholesky" )
global const PCPBJACOBI = String( "pbjacobi" )
global const PCMAT = String( "mat" )
global const PCHYPRE = String( "hypre" )
global const PCPARMS = String( "parms" )
global const PCFIELDSPLIT = String( "fieldsplit" )
global const PCTFS = String( "tfs" )
global const PCML = String( "ml" )
global const PCGALERKIN = String( "galerkin" )
global const PCEXOTIC = String( "exotic" )
global const PCCP = String( "cp" )
global const PCBFBT = String( "bfbt" )
global const PCLSC = String( "lsc" )
global const PCPYTHON = String( "python" )
global const PCPFMG = String( "pfmg" )
global const PCSYSPFMG = String( "syspfmg" )
global const PCREDISTRIBUTE = String( "redistribute" )
global const PCSVD = String( "svd" )
global const PCGAMG = String( "gamg" )
global const PCSACUSP = String( "sacusp" )
global const PCSACUSPPOLY = String( "sacusppoly" )
global const PCBICGSTABCUSP = String( "bicgstabcusp" )
global const PCAINVCUSP = String( "ainvcusp" )
global const PCBDDC = String( "bddc" )
global const PCKACZMARZ = String( "kaczmarz" )

typealias MatAssemblyType Cint
global const MAT_FLUSH_ASSEMBLY = Cint( 1 )
global const MAT_FINAL_ASSEMBLY = Cint( 0 )

typealias VecType String
global const VECSEQ = String( "seq" )
global const VECMPI = String( "mpi" )
global const VECSTANDARD = String( "standard" )
global const VECSHARED = String( "shared" )
global const VECSEQCUSP = String( "seqcusp" )
global const VECMPICUSP = String( "mpicusp" )
global const VECCUSP = String( "cusp" )
global const VECSEQVIENNACL = String( "seqviennacl" )
global const VECMPIVIENNACL = String( "mpiviennacl" )
global const VECVIENNACL = String( "viennacl" )
global const VECSEQCUDA = String( "seqcuda" )
global const VECMPICUDA = String( "mpicuda" )
global const VECCUDA = String( "cuda" )
global const VECNEST = String( "nest" )

typealias PCSide Cint
global const PC_SIDE_DEFAULT = Cint( -1 )
global const PC_LEFT = Cint( 0 )
global const PC_RIGHT = Cint( 1 )
global const PC_SYMMETRIC = Cint( 2 )

typealias MatType String
global const MATSAME = String( "same" )
global const MATMAIJ = String( "maij" )
global const MATSEQMAIJ = String( "seqmaij" )
global const MATMPIMAIJ = String( "mpimaij" )
global const MATIS = String( "is" )
global const MATAIJ = String( "aij" )
global const MATSEQAIJ = String( "seqaij" )
global const MATMPIAIJ = String( "mpiaij" )
global const MATAIJCRL = String( "aijcrl" )
global const MATSEQAIJCRL = String( "seqaijcrl" )
global const MATMPIAIJCRL = String( "mpiaijcrl" )
global const MATAIJCUSP = String( "aijcusp" )
global const MATSEQAIJCUSP = String( "seqaijcusp" )
global const MATMPIAIJCUSP = String( "mpiaijcusp" )
global const MATAIJCUSPARSE = String( "aijcusparse" )
global const MATSEQAIJCUSPARSE = String( "seqaijcusparse" )
global const MATMPIAIJCUSPARSE = String( "mpiaijcusparse" )
global const MATAIJVIENNACL = String( "aijviennacl" )
global const MATSEQAIJVIENNACL = String( "seqaijviennacl" )
global const MATMPIAIJVIENNACL = String( "mpiaijviennacl" )
global const MATAIJPERM = String( "aijperm" )
global const MATSEQAIJPERM = String( "seqaijperm" )
global const MATMPIAIJPERM = String( "mpiaijperm" )
global const MATSHELL = String( "shell" )
global const MATDENSE = String( "dense" )
global const MATSEQDENSE = String( "seqdense" )
global const MATMPIDENSE = String( "mpidense" )
global const MATELEMENTAL = String( "elemental" )
global const MATBAIJ = String( "baij" )
global const MATSEQBAIJ = String( "seqbaij" )
global const MATMPIBAIJ = String( "mpibaij" )
global const MATMPIADJ = String( "mpiadj" )
global const MATSBAIJ = String( "sbaij" )
global const MATSEQSBAIJ = String( "seqsbaij" )
global const MATMPISBAIJ = String( "mpisbaij" )
global const MATSEQBSTRM = String( "seqbstrm" )
global const MATMPIBSTRM = String( "mpibstrm" )
global const MATBSTRM = String( "bstrm" )
global const MATSEQSBSTRM = String( "seqsbstrm" )
global const MATMPISBSTRM = String( "mpisbstrm" )
global const MATSBSTRM = String( "sbstrm" )
global const MATDAAD = String( "daad" )
global const MATMFFD = String( "mffd" )
global const MATNORMAL = String( "normal" )
global const MATLRC = String( "lrc" )
global const MATSCATTER = String( "scatter" )
global const MATBLOCKMAT = String( "blockmat" )
global const MATCOMPOSITE = String( "composite" )
global const MATFFT = String( "fft" )
global const MATFFTW = String( "fftw" )
global const MATSEQCUFFT = String( "seqcufft" )
global const MATTRANSPOSEMAT = String( "transpose" )
global const MATSCHURCOMPLEMENT = String( "schurcomplement" )
global const MATPYTHON = String( "python" )
global const MATHYPRESTRUCT = String( "hyprestruct" )
global const MATHYPRESSTRUCT = String( "hypresstruct" )
global const MATSUBMATRIX = String( "submatrix" )
global const MATLOCALREF = String( "localref" )
global const MATNEST = String( "nest" )

typealias NormType Int32
global const NORM_1 = Int32( 0 )
global const NORM_2 = Int32( 1 )
global const NORM_FROBENIUS = Int32( 2 )
global const NORM_INFINITY = Int32( 3 )
global const NORM_MAX = Int32( 3 )

typealias PetscMatStructure Int32
global const DIFFERENT_NONZERO_PATTERN = Int32( 0 )
global const SUBSET_NONZERO_PATTERN = Int32( 1 )
global const SAME_NONZERO_PATTERN = Int32( 2 )

typealias MatReuse Int32
global const MAT_INITIAL_MATRIX = Int32( 0 )
global const MAT_REUSE_MATRIX = Int32( 1 )
global const MAT_IGNORE_MATRIX = Int32( 2 )

