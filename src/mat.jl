export PetscMat, MatSetType, SetUp, MatSetValues, MatAssemblyBegin, MatAssemblyEnd, MatSetSizes, MatGetSize, MatGetValues, MatGetOwnershipRange, MatXAIJSetPreallocation, MatMPIAIJSetPreallocation, MatSetFromOptions, MatGetInfo, MatMatMult, MatNorm, MatZeroEntries, MatSetValuesBlocked, MatSetOption, MatCreateShell, MatShellSetOperation, MatShellGetContext, MatGetType, MatCreateTranspose, MatTranspose


type PetscMat  <: AbstractArray{PetscScalar, 2}
  pobj::Ptr{Void}
  function PetscMat(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall( (:MatCreate,  libpetsclocation), PetscErrorCode, (comm_type, Ptr{Ptr{Void}}),comm,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection?
    return vec
  end

  function PetscMat(pobj::Ptr{Void})  # default constructor
    return new(pobj)
  end
end

#=
function PetscMat(mglobal::Integer, nglobal::Integer, comm::MPI_comm; mlocal=PETSC_DECIDE, nlocal=PETSC_DECIDE)

  mat = PetscMat(comm)
  PetscMatSetSize(
=#

function MatCreateShell(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer, arg6::Ptr{Void})
  # arg6 is the user provided context
    arg7 = Ref{Ptr{Void}}()
    ccall((:MatCreateShell,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Void}, Ref{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return PetscMat(arg7[])
end

function MatShellSetOperation(arg1::PetscMat,arg2::MatOperation,arg3::Ptr{Void})
# arg3 is a function pointer, and must have the signature:
# void fname(Mat, vec, vec) for MATOP_MULT

    ccall((:MatShellSetOperation,petsc),PetscErrorCode,(Ptr{Void},MatOperation,Ptr{Void}),arg1.pobj,arg2,arg3)
end

function MatShellGetContext(arg1::PetscMat)
# get the user provided context for the matrix shell
    arg2 = Ref{Ptr{Void}}()
    ccall((:MatShellGetContext,petsc),PetscErrorCode,(Ptr{Void},Ref{Ptr{Void}}),arg1.pobj,arg2)
    return arg2[]  # turn it into a julia object here?
end

function MatGetType(arg1::PetscMat)
    arg2 = Ref{Ptr{UInt8}}()
    ccall((:MatGetType,petsc),PetscErrorCode,(Ptr{Void}, Ref{Ptr{UInt8}}),arg1.pobj,arg2)
    return bytestring(arg2[])
end

function MatCreateTranspose(arg1::PetscMat)
  arg2 = Ref{Ptr{Void}}()

  ccall((:MatCreateTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Ptr{Void}}),arg1.pobj, arg2)

  return PetscMat(arg2[])
end

"""
  Constructs the transpose of a given matrix.  Can be in place or out of place.

  Inputs:
    A: matrix to take the transpose of

  Keywords:
    inplace: whether or not to do the transpose in place

  Outputs:
    A PetscMat object, either the original object A if the transpose was done 
    in place, or a new matrix object if the transpose was done out of place
"""
function MatTranspose(A::PetscMat; inplace::Bool=false, mat_initialized=false)
  if inplace
    B = Ref{Ptr{Void}}(A.pobj)
    reuse = MAT_REUSE_MATRIX
  else
    B = Ref{Ptr{Void}}()
    reuse = MAT_INITIAL_MATRIX
  end

  println("before, A.pobj = ", A.pobj)
    ccall((:MatTranspose,petsc),PetscErrorCode,(Ptr{Void},MatReuse,Ptr{Ptr{Void}}),A.pobj, reuse, B)

  println("after, A.pobj = ", A.pobj)
  println("typeof(B) = ", typeof(B))
  println("after, B[] = ", B[])

  if inplace
    return A
  else
    return PetscMat(B[])
  end
end






  function PetscDestroy(vec::PetscMat)
    if (vec.pobj != 0)
      err = ccall( (:MatDestroy,  libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
    end
#    vec.pobj = 0

#    println("Petsc Mat Destroy called")
#    sleep(5)
  end


immutable MatInfo
    block_size::PetscLogDouble
    nz_allocated::PetscLogDouble
    nz_used::PetscLogDouble
    nz_unneeded::PetscLogDouble
    memory::PetscLogDouble
    assemblies::PetscLogDouble
    mallocs::PetscLogDouble
    fill_ratio_given::PetscLogDouble
    fill_ratio_needed::PetscLogDouble
    factor_mallocs::PetscLogDouble

    function MatInfo()  # incomplete initialization
      return new()
    end
end


function show(io::IO, obj::MatInfo)
# print the fields of PetscMatInfo
#  println("PetscMatInfo:")
  println(io, "  block_size : ", obj.block_size)
  println(io, "  nz_allocated : ", obj.nz_allocated)
  println(io, "  nz_used : ", obj.nz_used)
  println(io, "  nz_unneeded : ", obj.nz_unneeded)
  println(io, "  memory : ", obj.memory)
  println(io, "  assemblies : ", obj.assemblies)
  println(io, "  mallocs : ", obj.mallocs)
  println(io, "  fill_ratio_given : ", obj.fill_ratio_given)
  println(io, "  fill_ratio_needed : ", obj.fill_ratio_needed)
  println(io, "  factor_mallocs : ", obj.factor_mallocs)
end

  function MatSetFromOptions(mat::PetscMat)
    ccall((:MatSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},), mat.pobj)

  end


  function MatSetType(vec::PetscMat,name)
    err = ccall( (:MatSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void}, Cstring), vec.pobj,name);

  end

  function SetUp(vec::PetscMat)
    err = ccall( ( :MatSetUp,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);

    MatSetOption(vec, MAT_ROW_ORIENTED, PETSC_TRUE)  # julia data is column-major
  end

#=
  PETSC_MAT_FLUSH_ASSEMBLY = 1;
  PETSC_MAT_FINAL_ASSEMBLY = 0
=#

  function MatSetValues(vec::PetscMat,idi::Array{PetscInt},idj::Array{PetscInt},array::Array{PetscScalar},flag::Integer)
    # remember, only matrices can be inserted into a Petsc matrix
    # if array is a 3 by 3, then idi and idj are vectors of length 3
#    idi = idi
#    idj = idj

    @assert length(idi)*length(idj) == length(array)

    # do check here to ensure array is the right shape (remember tranpose)
    err = ccall( ( :MatSetValues,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, Ptr{PetscInt}, PetscInt, Ptr{PetscInt}, Ptr{PetscScalar},Int32), vec.pobj,length(idi), idi, length(idj), idj,array,flag);
#    idi = idi
#    idj = idj
    return err
  end

  function MatSetValuesBlocked(mat::PetscMat, idi::Array{PetscInt}, idj::Array{PetscInt}, v::Array{PetscScalar}, flag::Integer)

      err = ccall((:MatSetValuesBlocked, libpetsclocation), PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}, Int32),mat.pobj, length(idi), idi, length(idj), idj, v,flag)

      return err
  end

  function MatSetOption(mat::PetscMat,arg2::MatOption,arg3::Union{Integer, Bool})
      ierr = ccall((:MatSetOption, libpetsclocation),PetscErrorCode, (Ptr{Void},MatOption,PetscBool), mat.pobj, arg2, arg3)

      if ierr != 0
	println(STDERR, "Error: MatSetOption returned non zero exit status")
      end
  end

  # 1-based indexing for both Petsc and regular matrices
  #----------------------------------------------------------------------------
  """
    1-based indexing for both regular and Pets matrices.
    The flag must be either PETSC_INSERT_VALUES or PETSC_ADD_VALUES.

    idxm and idxn cannot alias
  """
  function set_values1!(mat::PetscMat, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                        vals::Array{PetscScalar}, flag::Integer=PETSC_INSERT_VALUES)

    for i=1:length(idxm)
      idxm[i] -= 1
    end

    for i=1:length(idxn)
      idxn[i] -= 1
    end

    err = MatSetValues(mat, idxm, idxn, vals, flag)

    for i=1:length(idxm)
      idxm[i] += 1
    end

    for i=1:length(idxn)
      idxn[i] += 1
    end

    return err
  end

  function set_values1!{T}(mat::AbstractMatrix, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                           vals::Array{T}, flag::Integer=PETSC_INSERT_VALUES)

    if flag == PETSC_INSERT_VALUES
      for i=1:length(idxn)
        for j=1:length(idxm)
          mat[idxm[j], idxn[i]] = vals[j, i]
        end
      end
    elseif flag == PETSC_ADD_VALUES
      for i=1:length(idxn)
        for j=1:length(idxm)
          mat[idxm[j], idxn[i]] += vals[j, i]
        end
      end
    end

    return PetscErrorCode(0)
  end

  function set_values1!{T}(mat::SparseMatrixCSC, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                           vals::Array{T}, flag::Integer=PETSC_INSERT_VALUES)

    if flag == PETSC_INSERT_VALUES
      for i=1:length(idxn)
        for j=1:length(idxm)
          mat[idxm[j], idxn[i]] = vals[j, i]
        end
      end
    elseif flag == PETSC_ADD_VALUES  # optimized += implementation
      # hoist
      colptr = mat.colptr
      rowval = mat.rowval
      nzval = mat.nzval
      for i=1:length(idxn)
        row_start = colptr[idxn[i]]
        row_end = colptr[idxn[i]+1] - 1
        rowvals_extract = unsafe_view(rowval, row_start:row_end)
        for j=1:length(idxm)
          idx = searchsortedfirst(rowvals_extract, idxm[j])
          idx = row_start + idx - 1
          nzval[idx] += vals[j, i]
        end
      end
    end

    return PetscErrorCode(0)
  end


  function get_values1!(mat::PetscMat, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                        vals::Array{PetscScalar})

    for i=1:length(idxm)
      idxm[i] -= 1
    end

    for i=1:length(idxn)
      idxn[i] -= 1
    end

    err = MatGetValues(mat, idxm, idxn, vals)

    for i=1:length(idxm)
      idxm[i] += 1
    end

    for i=1:length(idxn)
      idxn[i] += 1
    end

    return err
  end

  function get_values1!(mat::AbstractMatrix, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                        vals::Array)

    for i=1:length(idxn)
      for j=1:length(idxm)
        vals[j, i] = mat[idxm[j], idxn[i]]
      end
    end

    return PetscErrorCode(0)
  end






  function MatAssemblyBegin(obj::AbstractArray, flg::Integer=0)
  end

  function MatAssemblyBegin(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
  end

  function MatAssemblyBegin(obj::PetscMat)
    return MatAssemblyBegin(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function MatAssemblyEnd(obj::AbstractArray, flag=0)
  end

  function MatAssemblyEnd(obj::PetscMat,flg::Integer)
    err = ccall( ( :MatAssemblyEnd,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
  end

  function MatAssemblyEnd(obj::PetscMat)
    return MatAssemblyEnd(obj,PETSC_MAT_FINAL_ASSEMBLY);
  end

  function MatSetSizes(vec::PetscMat,m::Integer, n::Integer, M::Integer, N::Integer)
    err = ccall( ( :MatSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void},PetscInt, PetscInt, PetscInt, PetscInt), vec.pobj,m,n,M,N);


    return err
  end

  function PetscView(obj::PetscMat,viewer)
    err = ccall( (:MatView,  libpetsclocation), PetscErrorCode, (Ptr{Void}, Int64),obj.pobj,0);
  end

  function MatGetSize(obj::PetscMat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    err = ccall(Libdl.dlsym(libpetsc, :MatGetSize), PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscInt}), obj.pobj,m,n);
    return m[1],n[1]
  end


export MatGetLocalSize
function MatGetLocalSize(mat::PetscMat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    ccall((:MatGetLocalSize,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, m, n)
    return m[1], n[1]
end




# Petsc populates v row major, so Julia will see the returned array as being transposed
function MatGetValues(obj::PetscMat, idxm::Array{PetscInt, 1}, idxn::Array{PetscInt, 1}, v::Array{PetscScalar, 2})
    # do check here to ensure v is the right shape

    ccall((:MatGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}), obj.pobj, length(idxm), idxm, length(idxn), idxn, v)

end

function MatGetOwnershipRange(mat::PetscMat)
    low = Array(PetscInt,1)
    high = Array(PetscInt,1)
    ccall((:MatGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, low, high)

    return low[1], high[1]
end


function MatZeroEntries(mat::PetscMat)
    ccall((:MatZeroEntries,petsc),PetscErrorCode,(Ptr{Void},),mat.pobj)
end

function MatZeroEntries(mat::AbstractMatrix)
  fill!(mat, 0.0)
end

function MatZeroEntries(mat::SparseMatrixCSC)
  fill!(mat.nzval, 0.0)
end



### new function ###

function MatNorm(mat::PetscMat, ntype::NormType)
    nrm = Array(PetscReal, 1)
    ccall((:MatNorm,petsc),PetscErrorCode,(Ptr{Void}, NormType,Ptr{PetscReal}), mat.pobj, ntype, nrm)
    return nrm[1]
end



export MatAXPY, MatAYPX, MatScale, MatShift
function MatAXPY(Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},PetscMatStructure), Y.pobj, a, X.pobj, str)
end

function MatAYPX( Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}, PetscMatStructure), Y.pobj, a, X.pobj, str)
end

function MatScale(mat::PetscMat, a::Number)
    ccall((:MatScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, PetscScalar(a))
end

function MatScale(mat::AbstractArray, a::Number)
  scale!(mat, PetscScalar(a))
end

function MatScale(mat::SparseMatrixCSC, a::Number)
  scale!(mat.nzval, PetscScalar(a))
end

function MatShift(mat::PetscMat, a::PetscScalar)
    ccall((:MatShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, a)
end

export MatMult, MatMultAdd, MatMultTranspose, MatMultHermitianTranspose

function MatMult(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMult,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), mat.pobj, x.pobj, y.pobj)
end


function MatMultAdd(A::PetscMat, x::PetscVec, y::PetscVec, z::PetscVec)
    ccall((:MatMultAdd,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj, z.pobj)
end

function MatMultTranspose(A::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj)
end

function MatMultHermitianTranspose(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultHermitianTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),mat.pobj, x.pobj, y.pobj)
end

function MatMatMult(A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat)
    arr = [C.pobj]
    ccall((:MatMatMult,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, MatReuse, PetscReal, Ptr{Ptr{Void}}), A.pobj, B.pobj, scall, fill, arr)
    if C.pobj != arr[1]
      PetscDestroy(C)
      C.pobj = arr[1] 
    end

end




function MatXAIJSetPreallocation(mat::PetscMat, bs::Integer, dnnz::AbstractArray{PetscInt, 1}, onnz::AbstractArray{PetscInt,1}, dnnzu::AbstractArray{PetscInt, 1}, onnzu::AbstractArray{PetscInt, 1})
# this is a unified interface for matrix preallocation for the Petsc built in
# matrix types: Aij, Bij, and their respective symmetric forms SAij, SBij
# for a non symmetric format matrix, dnnzu and onnzu are not required
# and vice versa for a non symmetric format matrix

    ccall((:MatXAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}), mat.pobj, bs, dnnz, onnz, dnnzu, onnzu)
end


# Matrix preallocation
function MatMPIAIJSetPreallocation(mat::PetscMat, d_nz::Integer, d_nnz::PetscInt_arr_or_null, o_nz::Integer, o_nnz::PetscInt_arr_or_null)

    ccall((:MatMPIAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void}, PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}), mat.pobj, d_nz, d_nnz, o_nz, o_nnz)
end

function MatGetInfo(mat::PetscMat, info_type::Int32)
    matinfo = MatInfo()  # create uninitialized struct
    ref_matinfo = Ref{MatInfo}(matinfo)
    ccall((:MatGetInfo,petsc),PetscErrorCode,(Ptr{Void}, Int32,Ref{MatInfo}), mat.pobj ,info_type, ref_matinfo)
    return ref_matinfo[]
end



# some functions to make a PetscMat look like a Julia array
# these only work after a matrix has been assembled
# some of them only work for uniprocess matrices
import Base.size
import Base.getindex

function size(A::PetscMat)

  return MatGetLocalSize(A)
end

function size(A::PetscMat, dim::Integer)
  dims = MatGetLocalSize(A)
  return dims[dim]
end

function getindex(A::PetscMat, i::Integer, j::Integer)
  # not efficient

  i_ = [PetscInt(i - 1)]
  j_ = [PetscInt(j - 1)]
  val = Array(PetscScalar, 1, 1)

  MatGetValues(A, i_, j_, val)

  return val[1]
end

function getindex(A::PetscMat, idx::Integer)
  # linear indexing
  m,n = MatGetLocalSize(A)
  i = idx % m
  j = div(idx, m)

  # zero based indexing
  i -= 1
  j -= 1  

  i_ = [PetscInt(i)]
  j_ = [PetscInt(j)]
  val = Array(PetscScalar, 1, 1)

  MatGetValues(A, i_, j_, val)
  return val[1]
end
