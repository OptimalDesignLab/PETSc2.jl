function test_options()
  @testset "----- testing options -----" begin

    PetscOptionsSetValue("-testkey", "42")
    PetscOptionsView()
    PetscOptionsClearValue("-testkey")
  end  # end facts
end 

test_options()
