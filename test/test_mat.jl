# test Petsc matrix functions
function make_mat()
   A = PetscMat(PETSC_DECIDE, PETSC_DECIDE,  "mpiaij", comm, mlocal=sys_size_local, nlocal=sys_size_local)
  SetUp(A)
  low, high = MatGetOwnershipRange(A)
  global_indices = Array(low:PetscInt(high-1))

  MatSetValues(A, global_indices, global_indices, A_julia, PETSC_INSERT_VALUES)
  MatAssemblyBegin(A)
  MatAssemblyEnd(A)

  return A
end

function make_vec()
  x = PetscVec(PETSC_DECIDE, VECMPI,  comm, mlocal=sys_size_local)

  low, high = VecGetOwnershipRange(x)
  global_indices = Array(low:PetscInt(high - 1))
 
  for i=1:sys_size_local
    idxm = [global_indices[i]]   # index
    val = [ rhs[i] ]  # value
    VecSetValues(x, idxm, val, PETSC_INSERT_VALUES)
  end

  VecAssemblyBegin(x)
  VecAssemblyEnd(x)

  return x
end



function test_indexing()
  facts("----- Testing Matrix indexing -----") do

    A = make_mat()
    B = make_mat()
    
    x = make_vec()
    low, high = VecGetOwnershipRange(x)
    global_indices = Array(low:PetscInt(high - 1))
    xvec_julia = copy(rhs)

    low, high = MatGetOwnershipRange(A)
    global_indices = Array(low:PetscInt(high - 1))
    
    B_julia = A_julia + PetscScalar(1)

    # test set_values1!
    global_indices1 = global_indices + PetscInt(1)
    global_indices2 = global_indices + PetscInt(1)
    vals = rand(PetscScalar, sys_size_local, sys_size_local)


    set_values1!(A, global_indices1, global_indices2, vals)
    MatAssemblyBegin(A)
    MatAssemblyEnd(A)
    vals2 = zeros(vals)
    MatGetValues(A, global_indices, global_indices, vals2)

    @fact vals --> roughly(vals2, atol=1e-13)
    fill!(vals2, 0.0)
    get_values1!(A, global_indices1, global_indices2, vals2)
    @fact vals --> roughly(vals2, atol=1e-13)


    A2s = zeros(3, 3)
    idxm = collect(PetscInt, 1:3)
    idxn = collect(PetscInt, 1:3)
    vals = rand(3,3)
    set_values1!(A2s, idxm, idxn, vals)

    @fact A2s[idxm,  idxn] --> roughly(vals, atol=1e-13)

    set_values1!(A2s, idxm, idxn, vals, PETSC_ADD_VALUES)

    @fact A2s[idxm, idxn] --> roughly(2*vals, atol=1e-13)
    vals2 = zeros(vals)
    get_values1!(A2s, idxm, idxn, vals2)
    @fact vals2 --> roughly(2*vals, atol=1e-13)


    
    # test the sparse optimized version of +=
    A2ss = sparse(A2s)
    set_values1!(A2ss, idxm, idxn, vals, PETSC_ADD_VALUES)
    fill!(vals2, 0.0)
    get_values1!(A2ss, idxm, idxn, vals2)
    @fact vals2 --> roughly(3*vals, atol=1e-13)

    for i=1:sys_size_local
      for j = 1:sys_size_local
        idxm = [ global_indices[i] ] # row index
        idxn = [ global_indices[j] ] # column index
        MatSetValues(A,idxm, idxn, [A_julia[i,j]],PETSC_INSERT_VALUES);
        MatSetValues(B,idxm, idxn, [A_julia[i,j] + PetscScalar(1)],PETSC_INSERT_VALUES);
      end
    end


    MatAssemblyBegin(A,PETSC_MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,PETSC_MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(B,PETSC_MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,PETSC_MAT_FINAL_ASSEMBLY);

    PetscDestroy(A)
    PetscDestroy(B)
    PetscDestroy(x)

  end  # end facts block

  return nothing
end  # end function


function test_transpose()
  facts("----- Testig Matrix Transpose -----") do

    A = make_mat()
    # test CreateMatTransposed
    At = MatCreateTranspose(A)

    b2t = make_vec()
    low, high = VecGetOwnershipRange(b2t)
    global_indices = Array(low:PetscInt(high - 1))
   
    x2 = make_vec()

    MatMult(At, x2, b2t)

    b2t_julia = A_julia*rhs

    vals = zeros(PetscScalar, sys_size_local)
    VecGetValues(b2t, global_indices, vals)

    for i=1:sys_size_local
      @fact b2t_julia[i] --> roughly(vals[i], atol=1e-12)
    end

    PetscDestroy(At)
    PetscDestroy(b2t)
    PetscDestroy(x2)


    # test MatTranpose
    At = MatTranspose(A, inplace=false)
    vals = zeros(PetscScalar, sys_size_local, sys_size_local)
    valst = zeros(vals)
    MatGetValues(A, global_indices, global_indices, vals)
    MatGetValues(At, global_indices, global_indices, valst)

    for i=1:sys_size_local
      for j=1:sys_size_local
        @fact vals[i, j] --> roughly(valst[j, i], atol=1e-12)
      end
    end

    PetscDestroy(At)

    At = MatTranspose(A, inplace=true)
  #  fill!(vals, 0.0)
    fill!(valst, 0.0)
  #  MatGetValues(A, global_indices, global_indices, vals)
    MatGetValues(At, global_indices, global_indices, valst)

    for i=1:sys_size_local
      for j=1:sys_size_local
        @fact vals[i, j] --> roughly(valst[j, i], atol=1e-12)
      end
    end

    # transpose A back to normal
    MatTranspose(A, inplace=true)
  #  PetscDestroy(At)

    PetscDestroy(A)
    PetscDestroy(b2t)
    PetscDestroy(x2)

  end

  return nothing
end


#=
  for i=1:sys_size_local
    for j=1:sys_size_local
      idxm = [global_indices[i] ]  # row index
      idxn = [ global_indices[j] ] # column index
      v = zeros(PetscScalar, 1,1)
      MatGetValues(A, idxm, idxn, v)
#      println("i = ", i, " j = ", j, " v = ", v[1,1], " A[i,j] = ", A_julia[i,j])
      @fact v[1,1] => roughly(A_julia[i,j]) "mismatch at i=$i, j=$j"
      @fact A[i,j] => roughly(A_julia[i,j])
    end
  end

  for i=1:sys_size_local*sys_size_local
    A[i] => roughly(A_julia[i])
  end
=#

function test_blas()
facts("----- Testing Matrix BLAS -----") do
  A = make_mat()
  B = make_mat()
  C = make_mat()

  # transpose everything to be column major
  MatTranspose(A, inplace=true)
  MatTranspose(B, inplace=true)
  MatTranspose(C, inplace=true)

  B_julia = copy(A_julia)
  B_copy = zeros(PetscScalar, sys_size_local, sys_size_local)

  y = make_vec()
  y_julia = copy(rhs)

  x = make_vec()
  x_julia = copy(rhs)
  xvec_julia = copy(rhs)

  low, high = VecGetOwnershipRange(x)
  global_indices = Array(low:PetscInt(high - 1))

  @fact size(A) => (sys_size_local, sys_size_local)
  @fact size(A, 1) => sys_size_local
  @fact size(A, 2) => sys_size_local


  @fact MatGetSize(A) => (comm_size*sys_size_local, comm_size*sys_size_local)
  @fact MatGetLocalSize(A) => (sys_size_local, sys_size_local)
  # testing non zero exist code is all we can really do here
  @fact PetscView(A) => 0

  alpha = PetscScalar(2.3)


  MatAXPY(B, alpha, A, DIFFERENT_NONZERO_PATTERN)
  B_julia = alpha*A_julia + B_julia
  MatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size_local
    for j=1:sys_size_local
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  MatAYPX(B, alpha, A, DIFFERENT_NONZERO_PATTERN)
  B_julia = alpha*B_julia + A_julia
  MatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size_local
    for j=1:sys_size_local
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  MatScale(B, alpha)
  B_julia = alpha*B_julia
  MatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size_local
    for j=1:sys_size_local
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  fac = 2.0
  C_julia = rand(3,3)
  C_julia2 = fac*C_julia
  scale!(C_julia, 2.0)
  @fact C_julia => roughly(C_julia2)

  C_julia = sprand(10, 10, 0.1)
  C_julia2 = fac*C_julia
  scale!(C_julia, fac)
  @fact C_julia => roughly(C_julia2)

  MatShift(B, alpha)
  B_julia += alpha*eye(PetscScalar, sys_size_local)
  MatGetValues(B, global_indices, global_indices, B_copy)
  for i=1:sys_size_local
    for j=1:sys_size_local
      @fact B_julia[j, i] => roughly(B_copy[i, j])
    end
  end

  y_copy = zeros(PetscScalar, sys_size_local)
  MatMult(B, x, y)
  y_julia = B_julia*xvec_julia

  VecGetValues(y, sys_size_local, global_indices, y_copy)
  for i=1:sys_size_local
      @fact y_julia[i] => roughly(y_copy[i])
  end


  MatMultAdd(B, x, y, y)
  y_julia = y_julia + B_julia*xvec_julia
  VecGetValues(y, sys_size_local, global_indices, y_copy)
  for i=1:sys_size_local
      @fact y_julia[i] => roughly(y_copy[i])
  end

  MatMultTranspose(A, x, y)
  y_julia = A_julia.'*xvec_julia
  VecGetValues(y, sys_size_local, global_indices, y_copy)
  for i=1:sys_size_local
      @fact y_julia[i] => roughly(y_copy[i])
  end

  MatMultHermitianTranspose(A, x, y)
  y_julia = A_julia'*xvec_julia
  VecGetValues(y, sys_size_local, global_indices, y_copy)
  for i=1:sys_size_local
      @fact y_julia[i] => roughly(y_copy[i])
  end


  D = PetscMat(comm)
  MatMatMult(A, B, MAT_INITIAL_MATRIX, PetscReal(1.0), D)
  D_julia = A_julia*B_julia
  D_copy = zeros(PetscScalar, sys_size_local, sys_size_local)
  MatGetValues(D, global_indices, global_indices, D_copy)

  for i=1:sys_size_local
    for j=1:sys_size_local
      @fact D_copy[j,i] => roughly(D_julia[i,j])
    end
  end

  # test matrix norms
  fnorm = MatNorm(D, NORM_FROBENIUS)
  infnorm = MatNorm(D, NORM_INFINITY)
  
  @fact fnorm => roughly(sqrt(comm_size)*vecnorm(D_julia), atol=1e-2)
  @fact infnorm => roughly(norm(D_julia, Inf))

  PetscDestroy(A)
  PetscDestroy(B)
  PetscDestroy(C)
  PetscDestroy(D)
  PetscDestroy(x)
  PetscDestroy(y)

end
return nothing
end



function test_prealloc()
facts("----- Testing Matrix Preallocation -----") do

  println("testing preallocation")
  nb = PetscInt(100)  # number of times/blocks to insert
  C = PetscMat(comm)
  MatSetType(C, "mpiaij")
#  MatSetFromOptions(C)

  MatSetSizes(C, nb*sys_size_local, nb*sys_size_local, PetscInt(nb*comm_size*sys_size_local),PetscInt(nb*comm_size*sys_size_local));

  # preallocation parameters
  bs = PetscInt(1)
  dnnz = 3*ones(PetscInt, nb*sys_size_local)  # on diagonal (row + column owned by this process)
  onnz = zeros(PetscInt, nb*sys_size_local)  # no off diagonal (column not owned by this process)
  dnnzu = Array(PetscInt, 0)  # this is not a symmetric matrix, so unused
  onnzu = Array(PetscInt, 0)  # this is not a symmetric matrix, so unused

  MatXAIJSetPreallocation(C, bs, dnnz, onnz, dnnzu, onnzu)

#  SetUp(C);


  A_julia_t = A_julia.'  # transpose because C is row major

  low, high = MatGetOwnershipRange(C)
  idi = zeros(PetscInt, sys_size_local)  # row indices
  idj = zeros(PetscInt, sys_size_local)  # column indices

  for i=1:nb

    # get indices
    pos = 1
    for j=1:sys_size_local # insert block on diagonal
	idi[pos] = low + sys_size_local*(i - 1) + j - 1
	idj[pos] = low + sys_size_local*(i - 1) + j - 1
	pos += 1
    end

    # add a random component
    for i=1:sys_size_local
      for j=1:sys_size_local
	A_julia_t[i,j] += rand()
      end
    end

#    println("idi = ", idi)
#    println("idj = ", idj)
    MatSetValues(C, idi, idj, A_julia_t, PETSC_INSERT_VALUES)
  end


  MatAssemblyBegin(C, PETSC_MAT_FINAL_ASSEMBLY)
  MatAssemblyEnd(C, PETSC_MAT_FINAL_ASSEMBLY)

  matinfo = MatGetInfo(C, MAT_LOCAL)

  @fact matinfo.mallocs => roughly(0.0)

#  PetscView(C, 0)

  println("finished testing preallocation")

  MatZeroEntries(C)
  cnorm = MatNorm(C, NORM_FROBENIUS)
  @fact cnorm => roughly(0.0, atol=1e-3)  # norm is zero iff matrix is zero

  @fact PetscDestroy(C) => 0

end
return nothing
end


test_indexing()
test_transpose()
test_blas()
test_prealloc()
