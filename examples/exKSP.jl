# Solve a simple 1D diffusion problem (without boundary conditions)
#
# ARGS[1] is the size of the system
#
# This script can be run with an arbitrary number of processes:
# julia ./exKSP.jl 40
# or
# mpirun -np 2 julia ./exKSP.jl 40
# or 
# mpirun -np 4 julia ./exKSP.jl 40
#
# and so on

using MPI
using PETSc2
PetscInitialize()

npts = parse(Int, ARGS[1])

myrank = MPI.Comm_rank(MPI.COMM_WORLD)
commsize = MPI.Comm_size(MPI.COMM_WORLD)

opts = Dict{ASCIIString, ASCIIString}(
  "-ksp_monitor" => "",
  "-malloc_dump" => "",
  "-pc_type" => "bjacobi",
  "-sub_pc_type" => "ilu",
  "-sub_pc_factor_levels" => "0",
  "-ksp_gmres_modifiedgramschmidt" => "",
)

PetscSetOptions(opts)

stencil = PetscScalar[1.0, -2.0, 1.0]
n_off_diag = div(length(stencil), 2)  # number of points off the main diagonal

# create matrix
A = PetscMat(npts, npts, MATMPIAIJ, MPI.COMM_WORLD)

# create the vectors and populate the right hand side
x = PetscVec(npts, VECMPI, MPI.COMM_WORLD)
b = PetscVec(npts, VECMPI, MPI.COMM_WORLD)


# preallocate A (this is a slight overestimate)
lcl_indices = local_indices(b)
npts_local = length(lcl_indices)

d_nnz = zeros(PetscInt, npts_local)
o_nnz = zeros(PetscInt, npts_local)

if commsize == 1
  fill!(d_nnz, length(stencil))
else
  o_nnz[1] = n_off_diag + 1
  fill!(d_nnz, length(stencil))
  o_nnz[end] = n_off_diag + 1
end

MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz)

# populate the matrix
idxm = PetscInt[0]
idxn = zeros(PetscInt, length(stencil))
if myrank == 0
  # do the first row
  idxm[1] = 1
  fill!(idxn, -1)  # negative entries are ignored
  start_idx = n_off_diag + 1
  idx = 1
  for i = start_idx:length(stencil)
    idxn[i] = idx
    idx += 1
  end

  set_values1!(A, idxm, idxn, stencil, PETSC_ADD_VALUES)
end

if myrank == (commsize-1)
  # do the last row
  idxm[1] = npts
  fill!(idxn, -1)

  idx = npts - (length(stencil) - n_off_diag) + 1
  for i=1:(length(stencil) - n_off_diag)
    idxn[i] = idx
    idx += 1
  end

  set_values1!(A, idxm, idxn, stencil, PETSC_ADD_VALUES)
end

# do everything else
first_row = lcl_indices[1]
last_row = lcl_indices[end]
if myrank == 0
  first_row += n_off_diag
end
if myrank == (commsize - 1)
  last_row -= n_off_diag
end

for i = first_row:last_row
  idxm[1] = i
  idxn[1] = i-1
  idxn[2] = i
  idxn[3] = i+1
  set_values1!(A, idxm, idxn, stencil, PETSC_ADD_VALUES)
end

# start assembly
assembly_begin(A, PETSC_MAT_FINAL_ASSEMBLY)

# populate right hand side
b_tmp = VecGetArray(b)
for i=1:length(b_tmp)
  b_tmp[i] = lcl_indices[i]  # arbitrary
end
VecRestoreArray(b, b_tmp)

# create the ksp solver
ksp = KSP(MPI.COMM_WORLD)
SetFromOptions(ksp)
SetOperators(ksp, A, A)
SetUp(ksp)
SetTolerances(ksp, 1e-12, 1e-12, PETSc2.PETSC_DEFAULT, PETSc2.PETSC_DEFAULT)


# end matrix assembly
assembly_end(A, PETSC_MAT_FINAL_ASSEMBLY)

# view the matrix
PetscView(A)

# solve the system
KSPSolve(ksp, b, x)

# view the solution
PetscView(x)

PetscDestroy(A)
PetscDestroy(x)
PetscDestroy(b)
PetscDestroy(ksp)
PetscFinalize()












