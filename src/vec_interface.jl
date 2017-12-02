# a support the Julia-ish interface when it maps cleanly onto a distributed
# vector

export set_values1!, get_values1!, assembly_begin, assembly_end, size_local,
       size_global, length_local, length_global, local_indices, fill_zero!,
       diagonal_shift!

"""
  Create a vector of a given size.  Users can specify either the global
  dimension or the local dimension
"""
function PetscVec(mglobal::Integer, comm::MPI_Comm; mlocal=PETSC_DECIDE)
  vec = PetscVec(comm)
  VecSetSizes(vec::PetscVec, mlocal, mglobal)

  return vec
end

"""
  Create a PetscVec, setting both the type and the format.  Users can specify
  either the local or global dimensions
"""
function PetscVec(mglobal::Integer, format, comm::MPI_Comm; mlocal=PETSC_DECIDE)
  vec = PetscVec(mglobal, comm, mlocal=mlocal)
  VecSetType(vec, format)
  return vec
end

"""
  PetscDestroy for AbstractVector.  No-op
"""
function PetscDestroy(vec::AbstractVector)
end

# 1-based indexing that unifies regular matrices and Petsc matrices
# TODO: generalize to AbstractArray when possible
#----------------------------------------------------------------------------

"""
  1-based indexing for both regular vectors and Petsc vector

  **Inputs**

   * vals: values to add/insert into the vector, must be length(idx)
   * flag: PETSC_INSERT_VALUES or PETSC_ADD_VALUES

  **Inputs/Outputs**

  * vec: the vector, can be a Petsc vector or a julia vector
  * idx: (global) indices to add/insert vals into
  
  idx is listed as input/output because it may be modified during the function.
  It will be returned to its original values when the function exits.
  This is necessary to accomodate Petscs zero-based indexing interface


"""
function set_values1!(vec::PetscVec, idx::Array{PetscInt}, vals::Array{PetscScalar}, flag::Integer=PETSC_INSERT_VALUES)
  for i=1:length(idx)
    idx[i] -= 1
  end

  err = VecSetValues(vec, idx, vals, flag)

  for i=1:length(idx)
    idx[i] += 1
  end

  return err
end

"""
  Method for AbstractVector
"""
function set_values1!{T}(vec::AbstractVector, idx::Array{PetscInt}, vals::Array{T}, flag::Integer=PETSC_INSERT_VALUES)
  if flag == PETSC_INSERT_VALUES
    for i in idx
      vec[i] = vals[i]
    end
  elseif flag == PETSC_ADD_VALUES
    for i in idx
      vec[i] += vals[i]
    end
  end

  return PetscErrorCode(0)
end


"""
  Like [`set_values1!`](@ref) but for retrieving values.  Note that Petsc
  only supports retrieving values from the local part of the vector

  **Inputs**

   * vec: a vector, can be a julia vector or a Petsc vector.
  

  **Inputs/Outputs**

   * idx: indices to retrieve
   * vals: array to put the values into

"""
function get_values1!(vec::PetscVec, idx::Array{PetscInt}, vals::Array{PetscScalar})
  for i=1:length(idx)
    idx[i] -= 1
  end

  err = VecGetValues(vec, idx, vals)

  for i=1:length(idx)
    idx[i] += 1
  end

  return err
end

"""
  Method for AbstractVector
"""
function get_values1!(vec::AbstractVector, idx::Array{PetscInt}, vals::Array)
  for i=1:length(idx)
    vals[i] = vec[idx[i]]
  end

  return PetscErrorCode(0)
end

"""
  Calls [`VecAssemblyBegin`](@ref).  No-op for Julia vectors

  **Inputs**
  
   * vec: AbstractVector
"""
function assembly_begin(vec::PetscVec)
  VecAssemblyBegin(vec)
end

function assembly_begin(vec::AbstractVector)
end

"""
  Calls [`VecAssemblyEnd`](@ref).  No-op for Julia vectors

  **Inputs**
  
   * vec: AbstractVector
"""
function assembly_end(vec::PetscVec)
  VecAssemblyEnd(vec)
end

function assembly_end(vec::AbstractVector)
end

"""
  Print non-Petsc vector to a given IO (a Julia IO, not a Petsc IO).  Defaults
  to printing to STDOUT.

  **Inputs**

   * b: AbstractVector
   * f: IO, defaults to STDOUT
"""
function PetscView(b::AbstractVector, f::IO=STDOUT)
  println(f, "b = \n", a)
end

"""
  Size of local part of vector

  **Inputs**

   * A: AbstractVector
"""
function size_local(A::PetscVec)
  return (VecGetLocalSize(A), )
end

function size_local(A::AbstractVector)
  return size(A)
end

"""
  Size of global vector
  
  **Inputs**

   * A: AbstractVector
"""
function size_global(A::PetscVec)
  return (VecGetSize(A), )
end

function size_global(A::AbstractVector)
  return size(A)
end

"""
  Length of local part of vector

  **Inputs**

   * A: AbstractVector
"""
function length_local(A::AllVectors)
  return size_local(A)[1]
end

"""
  Length of global vector

  **Inputs**

   * A: AbstractVector

"""
function length_global(A::AllVectors)
  return size_local(A)[1]
end

"""
  Returns a UnitRange containing the (1-based) global indices owned by this
  process.

  **Inputs**

   * A: AbstractVector

  **Outputs**

   * rng: UnitRange
"""
function local_indices(vec::PetscVec)

  low, high = VecGetOwnershipRange(vec)
  return Int(low+1):Int(high)
end

function local_indices(vec::AbstractVector)
  return 1:length(vec)  #TODO: update in future versions of Julia
end


"""
  Fill a vector with zeros

  **Inputs**

   * A: AbstractVector
"""
function fill_zero!(A::PetscVec)
  VecSet(A, PetscScalar(0.0))
end

function fill_zero!(A::AbstractVector)
  fill!(A, 0.0)
end

import Base.fill!

"""
  Base.fill! for Petsc vectors
"""
function fill!(A::PetscVec, a)
  VecSet(A, PetscScalar(a))
end

import Base.maximum

"""
  Base.maximum
"""
function maximum(x::PetscVec)
  mval, idx = VecMax(x)
  return mval
end

import Base.minimum

"""
  Base.minimum
"""
function minimum(x::PetscVec)
  mval, idx = VecMin(x)
  return mval
end

import Base.copy

"""
  Base.copy
"""
function copy(vec::PetscVec)
  b2 = VecDuplicate(vec)
  VecCopy(vec, b2)

  return b2
end

import Base.copy!
"""
  Base.copy!
"""
function copy!(dest::PetscVec, src::PetscVec)
  VecCopy(src, dest)
end


#TODO:
#      +, .*, -, (check for in-place versions),
#      


#TODO: extend blas routines

import Base.scale!
"""
  Base.scale!
"""
function scale!(vec::PetscVec, a::Number)
  _a = PetscScalar(a)
  VecScale(vec, _a)
end

"""
  Add a given value to all elements of the vector

  **Inputs**

   * A: AbstractVector
   * a: number to shift by
"""
function diagonal_shift!(A::PetscVec, a::Number)
  VecShift(A, PetscScalar(a))
end

function diagonal_shift!(A::AbstractVector, a::Number)
  @simd for i=1:length(A)
    A[i] += a
  end
end

import Base.norm
"""
  Base.norm
"""
function norm(A::PetscVec, p=2)
  if p == 1
    _p = NORM_1
  elseif p == 2
    _p = NORM_2
  elseif p == Inf
    _p = NORM_INFINITY
  else
    error("Petsc Vectors do not support norm p = $p")
  end

  VecNorm(A, _p)
end


import Base.dot
"""
  Base.dot

  Dot product where the *first* vector is conjugated.  This is is the reverse
  of VecDot, where the *second* vector is conjugated
"""
function dot(x::PetscVec, y::PetscVec)
  VecDot(y, x)
end

import Base.sum
"""
  Base.sum
"""
function sum(x::PetscVec)
  VecSum(x)
end


