using FactCheck
using PETSc2

@testset "  ---Checking Petsc data types---" begin

  @test ( PetscScalar )== Float64
  @test ( PetscReal )== Float64
  @test ( PetscInt )== Int64

end

FactCheck.exitstatus()
