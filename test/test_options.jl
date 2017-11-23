function test_options()
  facts("----- testing options -----") do

    PetscOptionsSetValue("-testkey", "42")
    PetscOptionsView()
    PetscOptionsClearValue("-testkey")
  end  # end facts
end 

test_options()
