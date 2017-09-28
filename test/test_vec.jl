function make_vec(format_j)

  b = PetscVec(sys_size_global, format_j, comm)
  low, high = VecGetOwnershipRange(b)
  global_indices = Array(low:PetscInt(high - 1))

  idxm = global_indices + PetscInt(1)  # 1 based indexing
  vals = Array(PetscScalar, sys_size_local)
  for i=1:sys_size_local
    vals[i] = rhs[i]
  end
  set_values1!(b, idxm, vals)
  PETSc.VecAssemble(b)

  return b
end

#------------------------------------------------------------------------------
# functions that run tests

# test indexing
function test_indexing(format_j)
  println("testing indexing")

  b = make_vec(format_j)
  #=
  b = PetscVec(sys_size_global, format_j, comm)
  idxm = global_indices + PetscInt(1)  # 1 based indexing
  vals = Array(PetscScalar, sys_size_local)
  for i=1:sys_size_local
    vals[i] = i
  end
  set_values1!(b, idxm, vals)
  PETSc.VecAssemble(b)
  =#
  vals2 = zeros(rhs)

  low, high = VecGetOwnershipRange(b)
  global_indices = Array(low:PetscInt(high - 1)) 
  idxm = global_indices + PetscInt(1)  # 1 based indexing

  VecGetValues(b, global_indices, vals2)

  @fact rhs --> roughly(vals2, atol=1e-13)
  fill!(vals2, 0.0)
  get_values1!(b,  idxm, vals2)
  @fact rhs --> roughly(vals2, atol=1e-13)

  b2s = zeros(3)
  idxm = collect(PetscInt, 1:3)
  vals = rand(3)
  set_values1!(b2s, idxm, vals)

  @fact b2s --> roughly(vals, atol=1e-13)

  set_values1!(b2s, idxm, vals, PETSC_ADD_VALUES)

  @fact b2s --> roughly(2*vals, atol=1e-13)
  vals2 = zeros(vals)
  get_values1!(b2s, idxm, vals2)

  @fact vals2 --> roughly(2*vals, atol=1e-13)




  for i=1:sys_size_local
    idxm = [global_indices[i]]   # index
    val = [ rhs[i] ]  # value
    VecSetValues(b, idxm, val, PETSC_INSERT_VALUES)
  end

  PETSc.VecAssemble(b)

  PetscDestroy(b)
  return nothing
end


# test copy, getArray
function test_copy(format_j)
  println("testing copy")
  b = make_vec(format_j)

  low, high = VecGetOwnershipRange(b)
  global_indices = Array(low:PetscInt(high - 1)) 

  # check that the vector was set/assembled correctly
  # check all the methods of copying/accesing a vector work
  b_copy = zeros(PetscScalar, sys_size_local)
  b2_copy = zeros(PetscScalar, sys_size_local)
#  idx = Array(0:2)  
#  idx = Array(PetscInt, 3)
#  for i=1:sys_size_local
#    idx[i] = global_indices[i]
#  end

  b_arr, ptr_arr = VecGetArray(b)
  b_arr_ro, ptr_arr2 = VecGetArrayRead(b)
  b2 = VecDuplicate(b)
  VecCopy(b, b2)
  
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  VecGetValues(b2, sys_size_local, global_indices, b2_copy)
  for i=1:sys_size_local
     @fact b_copy[i] => roughly(rhs[i])
     @fact b_arr[i] => roughly(rhs[i])
     @fact b_arr_ro[i] => roughly(rhs[i])
     @fact b2_copy[i] => roughly(rhs[i])
  end

  VecRestoreArray(b, ptr_arr)
  VecRestoreArrayRead(b, ptr_arr2)

 
  @fact VecGetSize(b) => comm_size*sys_size_local

  PetscDestroy(b)
  PetscDestroy(b2)
  return nothing
end

# test norms/BLAS stuff
function test_linalg(format_j)

  b = make_vec(format_j)
  b2 = copy(b)
#    b2 = VecDuplicate(b)
#    VecCopy(b, b2)
  # create 3rd vector to store results in
  b3 = PetscVec(sys_size_global, format_j, comm)

  low, high = VecGetOwnershipRange(b)
  global_indices = Array(low:PetscInt(high - 1)) 

  b_copy = zeros(PetscScalar, sys_size_local)
  b2_copy = zeros(PetscScalar, sys_size_local)
  rhs_tmp2 = zeros(PetscScalar, sys_size_local)
  for i=1:length(vec_norms)
    norm_petsc = VecNorm(b, vec_norms[i])
    norm_julia = norm(rhs, julia_vec_norms[i])
    @fact VecNorm(b, vec_norms[i]) => roughly( norm(rhs_global, julia_vec_norms[i]))
    print("\n")
  end
  
  # check for non zero exit status is all we can really do here 
  @fact PetscView(b, 0) => 0


  # test math functions
  VecSqrtAbs(b)
  rhs_tmp = deepcopy(rhs)  # don't modifiy original rhs
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = sqrt(abs(rhs_tmp[i]))
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecLog(b)
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = log(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecExp(b)
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = exp(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecAbs(b)
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = abs(rhs_tmp[i])
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  petsc_r, petsc_idx = VecMax(b)
  julia_r = maximum(real(b_copy))
  julia_idx = indmax(real(b_copy))

  @fact petsc_r => roughly(julia_r)
  @fact petsc_idx => (julia_idx - 1)

  petsc_r, petsc_idx = VecMin(b)
  julia_r = minimum(real(b_copy))
  julia_idx = indmin(real(b_copy))

  @fact petsc_r => roughly(julia_r)
  @fact petsc_idx => (julia_idx - 1)

  VecReciprocal(b)
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = 1/rhs_tmp[i]
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecShift(b, PetscScalar(2.0))
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp[i] = rhs_tmp[i] + 2.0
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecPointwiseMult(b3, b, b2)
  VecGetValues(b3, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp2[i] = rhs_tmp[i]*rhs[i]
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end


  VecPointwiseDivide(b3, b, b2)
  VecGetValues(b3, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    rhs_tmp2[i] = rhs_tmp[i]/rhs[i]
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end


  # test vector multiplication, addition functions
  # use b, b2 as the Petsc vectors
  # use rhs_tmp, rhs as the julia vectors

  alpha = PetscScalar(2.3)
  beta = PetscScalar(3.5)
  gamma = PetscScalar(4.8)

  VecAXPY(b, alpha, b2)
  rhs_tmp += alpha*rhs
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  #=
  VecAXPBY(b, alpha, beta, b2)
  rhs_tmp = alpha*rhs + beta*rhs_tmp
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end
=#

  VecAYPX(b, alpha, b2)
  rhs_tmp = alpha*rhs_tmp + rhs
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecWAXPY(b3, alpha, b, b2)
  rhs_tmp2 = alpha*rhs_tmp + rhs
  VecGetValues(b3, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp2[i])
  end

  # MAXPY
  scalar_arr = [alpha, beta]
  vec_arr = [b2.pobj, b3.pobj]
  VecMAXPY(b, PetscInt(2), scalar_arr, vec_arr)
  rhs_tmp += alpha*rhs + beta*rhs_tmp2
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecAXPBYPCZ(b, alpha, beta, gamma, b2, b3)
  rhs_tmp = alpha*rhs + beta*rhs_tmp2 + gamma*rhs_tmp
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end

  VecScale(b, alpha)
  rhs_tmp *= alpha
   VecGetValues(b, sys_size_local, global_indices, b_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs_tmp[i])
  end


  petsc_r = VecDot(b, b2)
  julia_r = comm_size*rhs'*rhs_tmp  # note reversed b, b2

  @fact petsc_r => roughly(julia_r[1])

#=
  petsc_r = VecTDot(b, b2)
  julia_r = rhs_tmp.'*rhs

  @fact petsc_r => roughly(julia_r[1])

  println("finished testing VecTDot")

  petsc_r =  VecSum(b2)
  julia_r = sum(rhs_global)

  @fact petsc_r => roughly(julia_r)
  
  println("finished testing VecSum")

  VecSwap(b, b2)
  b2_copy = zeros(PetscScalar, sys_size_local)
  VecGetValues(b, sys_size_local, global_indices, b_copy)
  VecGetValues(b2, sys_size_local, global_indices, b2_copy)
  for i=1:sys_size_local
    @fact b_copy[i] => roughly(rhs[i])
    @fact b2_copy[i] => roughly(rhs_tmp[i])
  end

  println("finished testing VecSwap")

  VecSwap(b, b2)  # swap back

=#

 
  
  PetscDestroy(b)
  PetscDestroy(b2)
  PetscDestroy(b3)
end




facts("\n   ---Testing vector functions---") do

for j=1:length(vec_formats)
#  format_j = "standard"
  format_j = vec_formats[j]
  println("testing vector format ", format_j)

  test_indexing(format_j)
  test_copy(format_j)
  test_linalg(format_j)

  end  # end loop over formats
end # end fact check block



