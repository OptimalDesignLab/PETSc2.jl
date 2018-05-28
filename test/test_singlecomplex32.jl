using FactCheck
using PETSc2
@testset "  ---Checking Petsc data types---" begin

  @test ( PetscScalar )== Complex64
  @test ( PetscReal )== Float32
  @test ( PetscInt )== Int32

end

FactCheck.exitstatus()
