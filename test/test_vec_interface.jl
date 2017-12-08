function test_vec_interface()
  # serial tests only
  if MPI.Comm_size(MPI.COMM_WORLD) != 1
    return nothing
  end
  facts("----- Testing Vec interface -----") do
    const m = 3
    vals = PetscScalar[1, 2, 3]
    vals2 = zeros(vals)  # used with get_values
    idx = PetscInt[1, 2, 3]
    x_j = zeros(PetscScalar, 3)
    set_values1!(x_j, idx, vals)

    assembly_begin(x_j)
    assembly_end(x_j)

    # test set_values!, both modes
    @fact norm(x_j - vals) --> roughly(0.0)
    set_values1!(x_j, idx, vals)
    @fact norm(x_j - vals) --> roughly(0.0)
    set_values1!(x_j, idx, vals, ADD_VALUES)
    @fact norm(x_j - 2*vals) --> roughly(0.0)
    get_values1!(x_j, idx, vals2)
    @fact norm(vals2 - 2*vals) --> roughly(0.0)

    x_p = PetscVec(3, PETSc2.VECMPI, MPI.COMM_WORLD)
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @fact norm(vals2 - vals) --> roughly(0.0)
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @fact norm(vals2 - vals) --> roughly(0.0)
    set_values1!(x_p, idx, vals, ADD_VALUES)
    assembly_begin(x_p)
    assembly_end(x_p)
    get_values1!(x_p, idx, vals2)
    @fact norm(vals2 - 2*vals) --> roughly(0.0)

    PetscView(x_p)  # if this doesn't error, that's all we can (easily) test

    @fact size_local(x_j) --> (3,)
    @fact size_global(x_j) --> (3,)
    @fact size_local(x_p) --> (3,)
    @fact size_global(x_p) --> (3,)
    @fact length_local(x_j) --> 3
    @fact length_global(x_j) --> 3
    @fact length_local(x_p) --> 3
    @fact length_global(x_p) --> 3
    @fact local_indices(x_j) --> 1:3
    @fact local_indices(x_p) --> 1:3


    fill_zero!(x_j)
    @fact norm(x_j) --> roughly(0.0)
    fill_zero!(x_p)
    @fact norm(x_p) --> roughly(0.0)

    fill!(x_j, 1.0)
    @fact norm(x_j, Inf) --> roughly(1.0)
    fill!(x_p, 1.0)
    @fact norm(x_p, Inf) --> roughly(1.0)

    set_values1!(x_j, idx, vals)
    set_values1!(x_p, idx, vals)
    assembly_begin(x_p)
    assembly_end(x_p)

    @fact maximum(x_p) --> roughly(3.0)
    @fact minimum(x_p) --> roughly(1.0)

    x2_p = copy(x_p)
    get_values1!(x2_p, idx, vals2)
    @fact norm(vals2 - vals) --> roughly(0.0)

    fill_zero!(x2_p)
    get_values1!(x2_p, idx, vals2)
    @fact norm(vals2) --> roughly(0.0)

    copy!(x2_p, x_p)
    get_values1!(x2_p, idx, vals2)
    @fact norm(vals2 - vals) --> roughly(0.0)



    PetscDestroy(x_p)
    PetscDestroy(x2_p)
  end
end

test_vec_interface()
