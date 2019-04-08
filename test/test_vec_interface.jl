function test_vec_interface()
  # serial tests only
  if MPI.Comm_size(MPI.COMM_WORLD) != 1
    return nothing
  end
  @testset "----- Testing Vec interface -----" begin
    const m = 3
    vals = PetscScalar[1, 2, 3]
    vals2 = zeros(vals)  # used with get_values
    idx = PetscInt[1, 2, 3]
    x_j = zeros(PetscScalar, 3)
    set_values1!(x_j, idx, vals)

    assembly_begin(x_j)
    assembly_end(x_j)

    # test set_values!, both modes
    @test isapprox( norm(x_j - vals), 0.0) 
    set_values1!(x_j, idx, vals)
    @test isapprox( norm(x_j - vals), 0.0) 
    set_values1!(x_j, idx, vals, ADD_VALUES)
    @test isapprox( norm(x_j - 2*vals), 0.0) 
    get_values1!(x_j, idx, vals2)
    @test isapprox( norm(vals2 - 2*vals), 0.0) 

    x_p = PetscVec(3, PETSc2.VECMPI, MPI.COMM_WORLD)
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @test isapprox( norm(vals2 - vals), 0.0) 
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @test isapprox( norm(vals2 - vals), 0.0) 
    set_values1!(x_p, idx, vals, ADD_VALUES)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @test isapprox( norm(vals2 - 2*vals), 0.0) 

    PetscView(x_p)  # if this doesn't error, that's all we can (easily) test

    @test ( size_local(x_j) )== (3,)
    @test ( size_global(x_j) )== (3,)
    @test ( size_local(x_p) )== (3,)
    @test ( size_global(x_p) )== (3,)
    @test ( length_local(x_j) )== 3
    @test ( length_global(x_j) )== 3
    @test ( length_local(x_p) )== 3
    @test ( length_global(x_p) )== 3
    @test ( local_indices(x_j) )== 1:3
    @test ( local_indices(x_p) )== 1:3


    fill_zero!(x_j)
    @test isapprox( norm(x_j), 0.0) 
    fill_zero!(x_p)
    @test isapprox( norm(x_p), 0.0) 

    fill!(x_j, 1.0)
    @test isapprox( norm(x_j, Inf), 1.0) 
    fill!(x_p, 1.0)
    @test isapprox( norm(x_p, Inf), 1.0) 

    set_values1!(x_j, idx, vals)
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)

    @test isapprox( maximum(x_p), 3.0) 
    @test isapprox( minimum(x_p), 1.0) 

    x2_p = copy(x_p)
    get_values1!(x2_p, idx, vals2)
    @test isapprox( norm(vals2 - vals), 0.0) 

    fill_zero!(x2_p)
    get_values1!(x2_p, idx, vals2)
    @test isapprox( norm(vals2), 0.0) 

    copy!(x2_p, x_p)
    get_values1!(x2_p, idx, vals2)
    @test isapprox( norm(vals2 - vals), 0.0) 

    # diagonal_shift
    x_j[:] = vals
    alpha = 2
    diagonal_shift!(x_j, alpha)
    @test isapprox( norm(x_j - (vals + alpha)), 0.0) atol=1e-14

    set_values1!(x_p, idx, vals, INSERT_VALUES)
    assembly_begin(x_p); assembly_end(x_p)
    diagonal_shift!(x_p, alpha)
    get_values1!(x_p, idx, vals2)
    @test isapprox( norm(vals2 - (vals + alpha)), 0.0) atol=1e-14

    PetscDestroy(x_p)
    PetscDestroy(x2_p)
  end
end

test_vec_interface()
