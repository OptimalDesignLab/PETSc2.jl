function test_mat_interface()

  if MPI.Comm_size(MPI.COMM_WORLD) != 1
    return nothing
  end

  @testset "----- Testing Mat interface -----" begin
    vals = PetscScalar[1 2 3; 4 5 6; 7 8 9]
    vals2 = zeros(vals)
    idx = PetscInt[1, 2, 3]
    idy = copy(idx)
    A_j = zeros(3,3)
    A_s = sparse(A_j)
    A_p = PetscMat(3, 3, PETSc2.MATMPIAIJ, MPI.COMM_WORLD)
    SetUp(A_p)

    test_setvalues1(A_j, idx, idy, vals)
    test_setvalues1(A_s, idx, idy, vals)
    test_setvalues1(A_p, idx, idy, vals)

    PetscView(A_j)
    PetscView(A_s)

    test_sizes(A_j)
    test_sizes(A_s)
    test_sizes(A_p)

    @test ( local_indices(A_j) )== 1:3
    @test ( local_indices(A_s) )== 1:3
    @test ( local_indices(A_p) )== 1:3

    fill_zero!(A_j)
    @test isapprox( norm(A_j, Inf), 0.0) 
    fill_zero!(A_s)
    @test isapprox( norm(A_s, Inf), 0.0) 
    fill_zero!(A_p)
    @test isapprox( norm(A_p, Inf), 0.0) 

    set_values1!(A_p, idx, idy, vals)
    assembly_begin(A_p, MAT_FINAL_ASSEMBLY)
    assembly_end(A_p, MAT_FINAL_ASSEMBLY)
    scale!(A_p, 2.0)
    get_values1!(A_p, idx, idy, vals2)
    @test isapprox( norm(vals2 - 2*vals), 0.0) atol=1e-14

    set_values1!(A_p, idx, idy, vals)
    assembly_begin(A_p, MAT_FINAL_ASSEMBLY)
    assembly_end(A_p, MAT_FINAL_ASSEMBLY)

    @test isapprox( vecnorm(A_p), vecnorm(vals)) atol=1e-13

    # test multiplication
    set_values1!(A_p, idx, idy, vals.')
    assembly_begin(A_p, MAT_FINAL_ASSEMBLY)
    assembly_end(A_p, MAT_FINAL_ASSEMBLY)


    x = PetscVec(3, VECMPI, MPI.COMM_WORLD)
    b = PetscVec(3, VECMPI, MPI.COMM_WORLD)
    vals_vec = PetscScalar[1, 2, 3]
    vals_vec2 = zeros(vals_vec)
    set_values1!(x, idx, vals_vec)

    A_mul_B!(b, A_p, x)
    get_values1!(b, idx, vals_vec2)
    b2 = vals*vals_vec
    @test isapprox( norm(vals_vec2 - b2), 0.0) atol=1e-14

    At_mul_B!(b, A_p, x)
    get_values1!(b, idx, vals_vec2)
    b2 = vals.'*vals_vec

    @test isapprox( norm(vals_vec2 - b2), 0.0) atol=1e-14

    PetscDestroy(A_p)
    PetscDestroy(x)
    PetscDestroy(b)
  end

  return nothing
end

function test_setvalues1(A, idx, idy, vals)

  vals2 = zeros(vals)

  # check set_values1!, both modes
  set_values1!(A, idx, idy, vals)
  assembly_begin(A, MAT_FINAL_ASSEMBLY)
  assembly_end(A, MAT_FINAL_ASSEMBLY)
  get_values1!(A, idx, idy, vals2)
  @test isapprox( norm(vals - vals2), 0.0) 

  set_values1!(A, idx, idy, vals)
  assembly_begin(A, MAT_FINAL_ASSEMBLY)
  assembly_end(A, MAT_FINAL_ASSEMBLY)
  get_values1!(A, idx, idy, vals2)
  @test isapprox( norm(vals - vals2), 0.0) 

  set_values1!(A, idx, idy, vals, ADD_VALUES)
  assembly_begin(A, MAT_FINAL_ASSEMBLY)
  assembly_end(A, MAT_FINAL_ASSEMBLY)
  get_values1!(A, idx, idy, vals2)
  @test isapprox( norm(2*vals - vals2), 0.0) 

  return nothing
end

function test_sizes(A)

  @test ( size_local(A) )== (3, 3)
  @test ( size_global(A) )== (3, 3)

  return nothing
end

test_mat_interface()
