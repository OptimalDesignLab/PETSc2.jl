export PetscMat, MatSetType, SetUp, MatSetValues, MatAssemblyBegin, MatAssemblyEnd, MatSetSizes, MatGetSize, MatGetValues, MatGetOwnershipRange, MatXAIJSetPreallocation, MatMPIAIJSetPreallocation, MatSetFromOptions, MatGetInfo, MatMatMult, MatNorm, MatZeroEntries, MatSetValuesBlocked, MatSetOption, MatCreateShell, MatShellSetOperation, MatShellGetContext, MatGetType, MatCreateTranspose, MatTranspose

"""
  PetscMat type.  Currently a subtype of `AbstractArray`, although that may
  change.
"""
type PetscMat  <: AbstractArray{PetscScalar, 2}
  pobj::Ptr{Void}

  """
    Constructor

    **Input**

     * comm: MPI communicator
  """
  function PetscMat(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall( (:MatCreate,  libpetsclocation), PetscErrorCode, (comm_type, Ptr{Ptr{Void}}),comm,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection?
    return vec
  end

  """
    Constructs a PetscMat object from an existing pointer to a Petsc matrix.

    **Inputs**

     * pobj: a void pointer to a Petsc matrix
  """
  function PetscMat(pobj::Ptr{Void})  # default constructor
    return new(pobj)
  end
end



"""
  MatCreateShell
"""
function MatCreateShell(arg1::MPI_Comm,arg2::Integer,arg3::Integer,arg4::Integer,arg5::Integer, arg6::Ptr{Void})
  # arg6 is the user provided context
    arg7 = Ref{Ptr{Void}}()
    ccall((:MatCreateShell,petsc),PetscErrorCode,(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Ptr{Void}, Ref{Ptr{Void}}),arg1,arg2,arg3,arg4,arg5,arg6,arg7)
    return PetscMat(arg7[])
end

"""
  MatShellSetOperation
"""
function MatShellSetOperation(arg1::PetscMat,arg2::MatOperation,arg3::Ptr{Void})
# arg3 is a function pointer, and must have the signature:
# void fname(Mat, vec, vec) for MATOP_MULT

    ccall((:MatShellSetOperation,petsc),PetscErrorCode,(Ptr{Void},MatOperation,Ptr{Void}),arg1.pobj,arg2,arg3)
end

"""
  MatShellGetContext
"""
function MatShellGetContext(arg1::PetscMat)
# get the user provided context for the matrix shell
    arg2 = Ref{Ptr{Void}}()
    ccall((:MatShellGetContext,petsc),PetscErrorCode,(Ptr{Void},Ref{Ptr{Void}}),arg1.pobj,arg2)
    return arg2[]  # turn it into a julia object here?
end




"""
  Frees a Petsc object.  Safe to call multiple times
"""
function PetscDestroy(vec::PetscMat)
  if (vec.pobj != 0)
    err = ccall( (:MatDestroy,  libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
  end
  vec.pobj = 0

#    println("Petsc Mat Destroy called")
#    sleep(5)
  return 0
end

"""
  Equivalent to Petsc's MatInfo struct
"""
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

import Base.show
"""
  show() for a MatInfo
"""
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

"""
  MatSetFromOptions
"""
function MatSetFromOptions(mat::PetscMat)
  ccall((:MatSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},), mat.pobj)

end

"""
  MatSetType
"""
function MatSetType(vec::PetscMat,name)
  err = ccall( (:MatSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void}, Cstring), vec.pobj,name);

end

"""
  MatSetUp
"""
function SetUp(vec::PetscMat)
  err = ccall( ( :MatSetUp,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);

end

"""
  MatGetType
"""
function MatGetType(arg1::PetscMat)
    arg2 = Ref{Ptr{UInt8}}()
    ccall((:MatGetType,petsc),PetscErrorCode,(Ptr{Void}, Ref{Ptr{UInt8}}),arg1.pobj,arg2)
    return bytestring(arg2[])
end



#=
  PETSC_MAT_FLUSH_ASSEMBLY = 1;
  PETSC_MAT_FINAL_ASSEMBLY = 0
=#

"""
  MatSetValues.  The `idx` and `idy` arrays must have PetscInt elements, and the
  `vals` array must have PetscScalar elements.
"""
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

"""
  MatSetValues blocked.
"""
function MatSetValuesBlocked(mat::PetscMat, idi::Array{PetscInt}, idj::Array{PetscInt}, v::Array{PetscScalar}, flag::Integer)

    err = ccall((:MatSetValuesBlocked, libpetsclocation), PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}, Int32),mat.pobj, length(idi), idi, length(idj), idj, v,flag)

    return err
end

"""
  MatSetOption
"""
function MatSetOption(mat::PetscMat,arg2::MatOption,arg3::Union{Integer, Bool})
    ierr = ccall((:MatSetOption, libpetsclocation),PetscErrorCode, (Ptr{Void},MatOption,PetscBool), mat.pobj, arg2, arg3)

    if ierr != 0
      println(STDERR, "Error: MatSetOption returned non zero exit status")
    end
end


"""
  MatAssemblyBegin
"""
function MatAssemblyBegin(obj::PetscMat,flg::Integer)
  err = ccall( ( :MatAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
end

"""
  MatAssemblyBegin, impicitly using `PETSC_MAT_FINAL_ASSEMBLY`
"""
function MatAssemblyBegin(obj::PetscMat)
  return MatAssemblyBegin(obj,PETSC_MAT_FINAL_ASSEMBLY);
end

"""
  MatAssemblyEnd
"""
function MatAssemblyEnd(obj::PetscMat,flg::Integer)
  err = ccall( ( :MatAssemblyEnd,  libpetsclocation), PetscErrorCode,(Ptr{Void},Int32), obj.pobj,flg);
end

"""
  MatAssemblyEnd, implicitly using `PETSC_MAT_FINAL_ASSEMBLY`
"""
function MatAssemblyEnd(obj::PetscMat)
  return MatAssemblyEnd(obj,PETSC_MAT_FINAL_ASSEMBLY);
end

"""
  MatSetSizes
"""
function MatSetSizes(vec::PetscMat,m::Integer, n::Integer, M::Integer, N::Integer)
  err = ccall( ( :MatSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void},PetscInt, PetscInt, PetscInt, PetscInt), vec.pobj,m,n,M,N);


  return err
end

"""
  Constructs the transpose of a given matrix.  Can be in place or out of place.

  **Inputs**

   * A: matrix to take the transpose of

  **Keywords**

   * inplace: whether or not to do the transpose in place

  **Outputs**
  
   * A PetscMat object, either the original object A if the transpose was done 
     in place, or a new matrix object if the transpose was done out of place
"""
function MatTranspose(A::PetscMat; inplace::Bool=false, mat_initialized=false)
  if inplace
    B = Ref{Ptr{Void}}(A.pobj)
    reuse = MAT_REUSE_MATRIX
  else
    B = Ref{Ptr{Void}}(C_NULL)
    reuse = MAT_INITIAL_MATRIX
  end

    ccall((:MatTranspose,petsc),PetscErrorCode,(Ptr{Void},MatReuse,Ptr{Ptr{Void}}),A.pobj, reuse, B)

  if inplace
    return A
  else
    return PetscMat(B[])
  end
end

"""
  MatCreateTranspose
"""
function MatCreateTranspose(arg1::PetscMat)
  arg2 = Ref{Ptr{Void}}()

  ccall((:MatCreateTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Ptr{Void}}),arg1.pobj, arg2)

  return PetscMat(arg2[])
end



"""
  PetscViewer
"""
function PetscView(obj::PetscMat,viewer::PetscViewer=C_NULL)
  err = ccall( (:MatView,  libpetsclocation), PetscErrorCode, (Ptr{Void}, Ptr{Void}),obj.pobj,viewer);
end

"""
  MatGetSize

  **Inputs**

   * obj: a PetscMat

  **Outputs**

   * m: first dimension size
   * n: second dimension size
"""
function MatGetSize(obj::PetscMat)
  m = Array(PetscInt, 1)
  n = Array(PetscInt, 1)
  err = ccall(Libdl.dlsym(libpetsc, :MatGetSize), PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscInt}), obj.pobj,m,n);
  return m[1],n[1]
end


export MatGetLocalSize

"""
  MatGetLocalsize.  Same interface as `MatGetSize`
"""
function MatGetLocalSize(mat::PetscMat)
    m = Array(PetscInt, 1)
    n = Array(PetscInt, 1)
    ccall((:MatGetLocalSize,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, m, n)
    return m[1], n[1]
end




# Petsc populates v row major, so Julia will see the returned array as being transposed
"""
  MatGetValues

  **Inputs**

   * obj: the PetscMat
   * idxm: row indices
   * idxn: column indices
   * vals: logically 2 dimensional array of values (although it could be a
           vector too)
"""
function MatGetValues(obj::PetscMat, idxm::Array{PetscInt, 1}, idxn::Array{PetscInt, 1}, v::Array{PetscScalar, 2})
    # do check here to ensure v is the right shape

    ccall((:MatGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}), obj.pobj, length(idxm), idxm, length(idxn), idxn, v)

end

"""
  MatGetOwnershipRange

  **Inputs**

   * mat: a PetscMat

  **Outputs**

   * low: lowest (zero-based) index that is owned
   * high: highest + 1 (zero-based) index that is owned
"""
function MatGetOwnershipRange(mat::PetscMat)
    low = Array(PetscInt,1)
    high = Array(PetscInt,1)
    ccall((:MatGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),mat.pobj, low, high)

    return low[1], high[1]
end

"""
  MatZeroEntries
"""
function MatZeroEntries(mat::PetscMat)
    ccall((:MatZeroEntries,petsc),PetscErrorCode,(Ptr{Void},),mat.pobj)
end



### new function ###

"""
  MatNorm
"""
function MatNorm(mat::PetscMat, ntype::NormType)
    nrm = Array(PetscReal, 1)
    ccall((:MatNorm,petsc),PetscErrorCode,(Ptr{Void}, NormType,Ptr{PetscReal}), mat.pobj, ntype, nrm)
    return nrm[1]
end



export MatAXPY, MatAYPX, MatScale, MatShift

"""
  MatAxPY
"""
function MatAXPY(Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},PetscMatStructure), Y.pobj, a, X.pobj, str)
end

"""
  MatAYPX
"""
function MatAYPX( Y::PetscMat, a::PetscScalar, X::PetscMat, str::PetscMatStructure)
    ccall((:MatAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}, PetscMatStructure), Y.pobj, a, X.pobj, str)
end

"""
  MatScale
"""
function MatScale(mat::PetscMat, a::Number)
    ccall((:MatScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, PetscScalar(a))
end

"""
  MatShift
"""
function MatShift(mat::PetscMat, a::Number)
    ccall((:MatShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), mat.pobj, a)
end

export MatMult, MatMultAdd, MatMultTranspose, MatMultHermitianTranspose

"""
  MatMult
"""
function MatMult(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMult,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), mat.pobj, x.pobj, y.pobj)
end

"""
  MatMultAdd
"""
function MatMultAdd(A::PetscMat, x::PetscVec, y::PetscVec, z::PetscVec)
    ccall((:MatMultAdd,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj, z.pobj)
end

"""
  MatMultTranspose
"""
function MatMultTranspose(A::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), A.pobj, x.pobj, y.pobj)
end

"""
  MatMultHermitianTranspose
"""
function MatMultHermitianTranspose(mat::PetscMat, x::PetscVec, y::PetscVec)
    ccall((:MatMultHermitianTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),mat.pobj, x.pobj, y.pobj)
end

"""
  MatMatMult
"""
function MatMatMult(A::PetscMat, B::PetscMat, scall::MatReuse, fill::PetscReal, C::PetscMat)
    arr = [C.pobj]
    ccall((:MatMatMult,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, MatReuse, PetscReal, Ptr{Ptr{Void}}), A.pobj, B.pobj, scall, fill, arr)
    if C.pobj != arr[1]
      PetscDestroy(C)
      C.pobj = arr[1] 
    end

end



"""
  MatXAIJSetPreallocation
"""
function MatXAIJSetPreallocation(mat::PetscMat, bs::Integer, dnnz::AbstractArray{PetscInt, 1}, onnz::AbstractArray{PetscInt,1}, dnnzu::AbstractArray{PetscInt, 1}, onnzu::AbstractArray{PetscInt, 1})
# this is a unified interface for matrix preallocation for the Petsc built in
# matrix types: Aij, Bij, and their respective symmetric forms SAij, SBij
# for a non symmetric format matrix, dnnzu and onnzu are not required
# and vice versa for a non symmetric format matrix

    ccall((:MatXAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt},Ptr{PetscInt}), mat.pobj, bs, dnnz, onnz, dnnzu, onnzu)
end


# Matrix preallocation
"""
  MatMPIAIJSetPreallocation
"""
function MatMPIAIJSetPreallocation(mat::PetscMat, d_nz::Integer, d_nnz::PetscInt_arr_or_null, o_nz::Integer, o_nnz::PetscInt_arr_or_null)

    ccall((:MatMPIAIJSetPreallocation,petsc),PetscErrorCode,(Ptr{Void}, PetscInt,Ptr{PetscInt},PetscInt,Ptr{PetscInt}), mat.pobj, d_nz, d_nnz, o_nz, o_nnz)
end

"""
  MatGetInfo
"""
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

include("mat_interface.jl")
