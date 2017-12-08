using MPI, PETSc2

  PetscInitialize()
  b = PetscVec(10, VECMPI, MPI.COMM_WORLD)

  idx = PetscInt[1, 2, 3]
  vals = PetscScalar[1.0, 2.0, 3.0]

  set_values1!(b, idx, vals, ADD_VALUES)

  # don't assemble the vector
  PetscView(b)

  PetscDestroy(b)
