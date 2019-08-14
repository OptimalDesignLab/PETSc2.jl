
using PETSc2
using Base.Test

function not_isapprox(args...; kwargs...)
  return !isapprox(args...; kwargs...)
end



import MPI
println("pwd = ", pwd())
println("PETSC_DIR = ", ENV["PETSC_DIR"])
if !MPI.Initialized()
  MPI.Init()
end

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)


# write your own tests here


# define some arrays to make looping easier
#vec_norms = [PETSC_NORM_1, PETSC_NORM_2, PETSC_NORM_INFINITY]

# skip 1 norm untill Petsc 1 complex 1 norm fix is released
vec_norms = [ NORM_2, NORM_INFINITY]
julia_vec_norms = [ 2, Inf]
# for some reason only the first four vector format works
#vec_formats = [VECSEQ, VECMPI, VECSTANDARD, VECSHARED, VECSEQCUSP, VECMPICUSP, VECCUSP, VECSEQVIENNACL, VECMPIVIENNACL, VECVIENNACL, VECNEST]  # no pthread

if comm_size > 1  # this is a parallel run
  vec_formats = [ VECMPI, VECSTANDARD]
else  # serial run 
  vec_formats = [VECSEQ]
end

PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug","-malloc_dump"]);
#PetscInitialize(["-ksp_monitor","-malloc","-malloc_debug"]);

# size of the system owned by this process (ie. local sys size)
sys_size_local = PetscInt(3)
sys_size_global = sys_size_local*comm_size

# indicies of the vector owned by this process
# create these with smallest precision, so they can be promoted
tmp2 = convert(Array{Float32,1}, collect(1.0:3))
tmp = convert(Array{Float32, 2}, [1.0 2.0 3; 4 5 7; 7 8 9])

tmp3 = convert(Array{Complex64, 1}, [1.0 + 0im; 2.0 + 1.0im; 3.0 + 2.0im])
tmp4 = convert( Array{Complex64, 2}, [1.0 + 1im   2 + 2im  3 + 3im; 4 + 4im  5 + 5im 7 + 7im; 7 + 7im 8 + 8im 9 + 9im])

# Have both local and global versions of system for testing
A_julia = zeros(PetscScalar, sys_size_local, sys_size_local)
rhs = zeros(PetscScalar, sys_size_local)
A_julia_global = zeros(PetscScalar, comm_size*sys_size_local, comm_size*sys_size_local)
rhs_global = zeros(PetscScalar, comm_size*sys_size_local)
# convert to arrays of the proper Petsc type
# this facilitates testing with the different Petsc build options

# define tolerances for KSP
if PetscReal == Float64
  rtol = PetscReal(1e-5)
  abstol = PetscReal(1e-14)
else  # single preciiosn
  rtol = PetscReal(1e-5)
  abstol = PetscReal(1e-6)
end
dtol = PetscReal(1e5)
maxits = PetscInt(50)


if PetscScalar <: Real

  for i=1:sys_size_local
    rhs[i] = convert(PetscScalar, tmp2[i])
    
    for j=0:(comm_size - 1)  # global 
      idx = i + j*sys_size_local
      rhs_global[idx] = convert(PetscScalar, tmp2[i])
    end    

  end

  for i=1:sys_size_local
    for j=1:sys_size_local
      A_julia[i,j] = convert(PetscScalar, tmp[i,j])

      for k=0:(comm_size - 1)
        idxm = i + k*sys_size_local
        idxn = j + k*sys_size_local
        A_julia_global[idxm, idxn] = convert(PetscScalar, tmp[i,j] )
      end
    end
  end

end  # end if statement

if PetscScalar <: Complex

  for i=1:sys_size_local
    rhs[i] = convert(PetscScalar, tmp3[i])

    for j=0:(comm_size - 1)  # global 
      idx = i + j*sys_size_local
      rhs_global[idx] = convert(PetscScalar, tmp3[i])
    end    


  end

  for i=1:sys_size_local
    for j=1:sys_size_local
      A_julia[i,j] = convert(PetscScalar, tmp4[i,j])

      for k=0:(comm_size - 1)
        idxm = i + k*sys_size_local
        idxn = j + k*sys_size_local
        A_julia_global[idxm, idxn] = convert(PetscScalar, tmp4[i,j] )
      end
 
    end
  end
end
#

if comm_rank == 0
  println("rhs_global = ", rhs_global)
  println("A_julia_global = ", A_julia_global)
end

#println("rhs = ", rhs)
#println("A_julia = ", A_julia)
x_julia = A_julia\rhs
#println("x_julia = ", x_julia)

include("test_vec.jl")
include("test_vec_interface.jl")
include("test_mat.jl")
include("test_mat_interface.jl")
include("test_ksp.jl")
include("test_options.jl")

#include("test_prealloc.jl")

PetscFinalize()
