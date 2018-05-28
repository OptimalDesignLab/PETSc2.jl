using FactCheck
using PETSc2
@testset "  ---Checking Petsc data types---" begin

  @test ( PetscScalar )== Complex64
  @test ( PetscReal )== Float32
  @test ( PetscInt )== Int64

end

FactCheck.exitstatus()
