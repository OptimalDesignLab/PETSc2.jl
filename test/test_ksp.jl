function test_ksp()

  @testset "\n   ---testing KSP solvers---" begin
  # create vectors and matrices

    b = make_vec()

    low, high = VecGetOwnershipRange(b)
    b_global_indices = collect(low:PetscInt(high - 1))
   
    x = make_vec()
    low, high = VecGetOwnershipRange(x)
    x_global_indices = collect(low:PetscInt(high - 1))

    println("making matrix for ksp solve")  
    A = make_mat()
    MatTranspose(A, inplace=true)

    low, high = MatGetOwnershipRange(A)
    mat_global_indices = collect(low:PetscInt(high - 1))
   

    # perform solve
    ksp = KSP(comm)
    SetOperators(ksp, A, A)
    SetFromOptions(ksp)
    SetUp(ksp)
    KSPSolve(ksp, b, x)
    reason = GetConvergedReason(ksp)

    @test  reason  > 0# convergence

    # copy solution back to Julia
    x_copy = zeros(PetscScalar, sys_size_local)
    VecGetValues(x, sys_size_local, x_global_indices, x_copy)
    for i=1:sys_size_local
        @test isapprox( x_copy[i], x_julia[i]) atol=1e-14
    end

    PetscDestroy(ksp)


    println("   \n--- Testing LGMRES --- ")
    # test using non default KSP method

    # perform solve
    ksp = KSP(comm)
    SetOperators(ksp, A, A)

    # test PC
    println("   \n---Testing PC Options ---")
    pc = KSPGetPC(ksp)

    PCSetType(pc, PETSc2.PCBJACOBI)
    pctype = PCGetType(pc)
    @test ( pctype )== PETSc2.PCBJACOBI
    #PCFactorSetUseInPlace(pc, PetscBool(true))

    ### Do some KSP setup

    SetFromOptions(ksp)
    SetTolerances(ksp, rtol, abstol, dtol, maxits)
    SetInitialGuessNonzero(ksp, PetscBool(true))
    SetType(ksp, PETSc2.KSPLGMRES)
    SetUp(ksp)

    ### More PC testing
    n_local, first_local, ksp_arr = PCBJacobiGetSubKSP(pc)
    ksp2 = ksp_arr[1]

    pc2 = KSPGetPC(ksp2)

    PCSetReusePreconditioner(pc2, PetscBool(true))
    @test ( PCGetReusePreconditioner(pc2) )== true

    PCFactorSetAllowDiagonalFill(pc2, PetscBool(true))
    @test ( PCFactorGetAllowDiagonalFill(pc2) )== true

    PCFactorSetLevels(pc2, PetscInt(1))
    @test ( PCFactorGetLevels(pc2) )== 1  # should be pc2

    PCSetReusePreconditioner(pc2, PetscBool(true))
    @test ( PCGetReusePreconditioner(pc2) )== true

    PCFactorSetFill(pc, PetscReal(7.0))

    #=
    PCJacobiSetType(pc, PETSc2.PC_JACOBI_ROWMAX)
    @test ( PCJacobiGetType(pc) )== PETSc2.PC_JACOBI_ROWMAX
    =#





    #=
    tmp = PCFactorGetUseInPlace(pc)
    println("tmp = ", tmp)
    @test ( PCFactorGetUseInPlace(pc) )== true
    =#


    KSPSolve(ksp, b, x)
    reason = GetConvergedReason(ksp)
    ksptype = GetType(ksp)

    @test ( ksptype )== PETSc2.KSPLGMRES


    rtol_ret, abstol_ret, dtol_ret, maxits_ret = GetTolerances(ksp)

    @test isapprox( rtol_ret, rtol) 
    @test isapprox( abstol_ret, abstol) 
    @test isapprox( dtol_ret, dtol) 
    @test ( maxits_ret )== maxits
    @test ( GetInitialGuessNonzero(ksp) )== true
    @test  reason  > 0# convergence

    rnorm = GetResidualNorm(ksp)
    @test  rnorm  < abstol

    # copy solution back to Julia
    x_copy = zeros(PetscScalar, sys_size_local)

    VecGetValues(x, sys_size_local, x_global_indices, x_copy)
    for i=1:sys_size_local
        @test isapprox( x_copy[i], x_julia[i]) atol=1e-14
    end





    PetscDestroy(x)
    PetscDestroy(b)
    PetscDestroy(A)
    PetscDestroy(ksp)
  end

return nothing
end

test_ksp()



