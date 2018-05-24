# interface for PetscMat that also maps cleanly onto the Julia Matrix interface
# the Julia sparse matrix interface is ill-defined, so this may take some
# liberties

#=
function PetscMat(mglobal::Integer, nglobal::Integer, comm::MPI_Comm; mlocal=PETSC_DECIDE, nlocal=PETSC_DECIDE)

  return mat
end
=#

"""
  Constructor

  **Inputs**

   * mglobal: first dimension global size (or PETSC_DECIDE)
   * nglobal: second dimension global size (or PETSC_DECIDE)
   * format: matrix format
   * comm: MPI communicator

  **Keyword Arguments**

   * mlocal: first dimension local size (or PETSC_DECIDE)
   * nlocal: second dimension local size (or PETSC_DECIDE)
"""
function PetscMat(mglobal::Integer, nglobal::Integer, format, comm::MPI_Comm;
                  mlocal=PETSC_DECIDE, nlocal=PETSC_DECIDE)

  mat = PetscMat(comm)
  MatSetSizes(mat, mlocal, nlocal, mglobal, nglobal)
  MatSetType(mat, format)


  return mat
end

#=
function PetscDestroy(vec::PetscMat)
  if (vec.pobj != 0)
    err = ccall( (:MatDestroy,  libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
  end
  vec.pobj = 0

#    println("Petsc Mat Destroy called")
#    sleep(5)
  return 0
end
=#

"""
  PetscDestroy for AbstractMatrix.  No-op
"""
function PetscDestroy(mat::AbstractMatrix)
end


# 1-based indexing for both Petsc and regular matrices
#----------------------------------------------------------------------------
"""
  1-based indexing for both regular and Pets matrices.
  Note that Petsc treats arrays at being row-major, so it is recommened
  to set MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE) before using
  this function.

  **Inputs**
  
   * flag: INSERT_VALUES or ADD_VALUES.  Note that the first one
           result in non-deterministic behavior in parallel (in the general
           case)

  * vals: the values, must be length(idxm) x length(idxn)

  **Inputs/Outputs**

   * mat: a matrix, can be a Petsc matrix or a julia matrix
   * idxm: the row numbers
   * idxn: the column numbers
   

  Note that idxm and idxn are listed as input/outputs because they may be
  modified by this function, however when the function returns they
  will have the same values as on entry.  This is needed to accomodate the
  fact that Petsc uses 1 based indexing internally.

  This function is optimized for PetscMat and SparseMatrixCSC

  Aliasing restriction: idxm and idxn cannot alias
"""
function set_values1!(mat::PetscMat, idxm::Array{PetscInt}, idxn::Array{PetscInt}, 
                      vals::Array{PetscScalar}, flag::Integer=INSERT_VALUES)

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

"""
  Method for AbstractMatrix
"""
function set_values1!{T}(mat::AbstractMatrix, idxm::Array{PetscInt},
                         idxn::Array{PetscInt}, vals::Array{T},
                         flag::Integer=INSERT_VALUES)

  if flag == INSERT_VALUES
    for i=1:length(idxn)
      for j=1:length(idxm)
        mat[idxm[j], idxn[i]] = vals[j, i]
      end
    end
  elseif flag == ADD_VALUES
    for i=1:length(idxn)
      for j=1:length(idxm)
        mat[idxm[j], idxn[i]] += vals[j, i]
      end
    end
  end

  return PetscErrorCode(0)
end

"""
  Method for `SparseMatrixCSC`
"""
function set_values1!{T}(mat::SparseMatrixCSC, idxm::Array{PetscInt},
                         idxn::Array{PetscInt}, vals::Array{T},
                         flag::Integer=INSERT_VALUES)

  if flag == INSERT_VALUES
    for i=1:length(idxn)
      for j=1:length(idxm)
        mat[idxm[j], idxn[i]] = vals[j, i]
      end
    end
  elseif flag == ADD_VALUES  # optimized += implementation
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

"""
  Like [`set_values1!`](@ref), but retrieves values.  See that function for
  the meanings of the arguments. Note that Petsc does
  not support getting values for the non-local block of the matrix

  **Inputs**

   * mat: a matrix, can be a Petsc matrix or a julia matrix

  **Inputs/Outputs**

   * idxm
   * idxn
   * vals

  Aliasing restrictions: idxm and idxn cannot alias
"""
function get_values1!(mat::PetscMat, idxm::Array{PetscInt},
                      idxn::Array{PetscInt}, 
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

"""
  Method for `AbstractMatrix`
"""
function get_values1!(mat::AbstractMatrix, idxm::Array{PetscInt},
                      idxn::Array{PetscInt}, vals::Array)

  for i=1:length(idxn)
    for j=1:length(idxm)
      vals[j, i] = mat[idxm[j], idxn[i]]
    end
  end

  return PetscErrorCode(0)
end

"""
  Begin matrix assembly for `PetscMat`.  No-op for Julia matrices

  **Inputs**

   * mat: AbstractMatrix
   * flag: type of matrix assembly (see [`MatAssemblyBegin`](@ref)
"""
function assembly_begin(mat::PetscMat, flg::Integer)
  MatAssemblyBegin(mat, flg)
end

function assembly_begin(mat::AbstractMatrix, flg::Integer)
end

"""
  Counterpart of [`assembly_end`](@ref)

  **Inputs**

   * mat: AbstractMatrix
   * flg: type of matrix assembly
"""
function assembly_end(mat::PetscMat, flg::Integer)
  MatAssemblyEnd(mat, flg)
end

function assembly_end(mat::AbstractMatrix, flg::Integer)
end

#TODO: transpose, inplace and/out out of place
import Base: transpose, transpose!

# PetscView?
"""
  Print a non-Petsc matrix to a given IO (a Julia IO, not a Petsc IO)

  **Inputs**

   * A: the matrix
   * f: an IO, default STDOUT
"""
function PetscView(A::AbstractMatrix, f::IO=STDOUT)
  println(f, "A = \n", A)
end

"""
  Size of local part of matrix

  **Inputs**

   * A: AbstractMatrix
"""
function size_local(A::PetscMat)
  return MatGetLocalSize(A)
end

function size_local(A::AbstractMatrix)
  return size(A)
end

"""
  Global size of matrix, same as size_local() for serial matrices

  **Inputs**

   * A: AbstractMatrix
"""
function size_global(A::PetscMat)
  return MatGetSize(A)
end

function size_global(A::AbstractMatrix)
  return size(A)
end

"""
  Returns the rows owned by this process (1-based)

  **Inputs**

   * A: a matrix

  **Outputs**

   * rng: a UnitRange
"""
function local_indices(A::PetscMat)
  low, high = MatGetOwnershipRange(A)
  return (low+1):high
end

function local_indices(A::AbstractMatrix)
  return 1:size(A, 1)
end

"""
  Fill the matrix with zeros. The sparsity pattern of the matrix (if applicable)
  should be defined before this function is called

  **Inputs**

   * A: AbstractMatrix
"""
function fill_zero!(A::PetscMat)
  MatZeroEntries(A)
end

function fill_zero!(A::AbstractMatrix)
  fill!(A, 0.0)
end

function fill_zero!(A::SparseMatrixCSC)
  fill!(A.nzval, 0.0)
end

import Base.scale!
"""
  scale! for Petsc matrix
"""
function scale!(A::PetscMat, a::Number)
  MatScale(A, PetscScalar(a))
end

"""
  Adds the specified value to the diagonal of the matrix

  **Inputs**

   * A: the matrix
   * a: the value
"""
function diagonal_shift!(A::PetscMat, a::Number)
  MatShift(A, PetscScalar(a))
end

function diagonal_shift!(A::AbstractMatrix, a::Number)

  nmax = min(size(A, 1), size(A, 2))
  @simd for i=1:nmax
    A[i, i] += 1
  end
end
#TODO; optimized implementation for SparseMatrixCSC


import Base: norm, vecnorm

"""
  Norm for Petsc matrices, 1, 2, and infinity norms supported
"""
function norm(A::PetscMat, p::Real)
  if p == 1
    _p = NORM_1
  elseif p == 2
    _p = NORM_2
  elseif p == Inf
    _p = NORM_INFINITY
  end
  MatNorm(A, _p)
end

"""
  Frobenius norm, consistent with Julias interface
"""
function vecnorm(A::PetscMat)
  MatNorm(A, NORM_FROBENIUS)
end

import Base: A_mul_B!, At_mul_B!

"""
  Computes b = A*x
"""
function A_mul_B!(b::PetscVec, A::PetscMat, x::PetscVec)
  MatMult(A, x, b)
  return b
end

"""
  Computes b = A.'*x
"""
function At_mul_B!(b::PetscVec, A::PetscMat, x::PetscVec)
  MatMultTranspose(A, x, b)
  return b
end

