# a support the Julia-ish interface when it maps cleanly onto a distributed
# vector

#TODO: set Petsc to interpret arrays as column-major

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


function get_values1!(vec::AbstractVector, idx::Array{PetscInt}, vals::Array)
  for i=1:length(idx)
    vals[i] = vec[idx[i]]
  end

  return PetscErrorCode(0)
end

#TODO: VecAssemblyBegin/End?
function assembly_begin(vec::PetscVec)
  VecAssemblyBegin(vec)
end

function assembly_begin(vec::AbstractVector)
end

function assembly_end(vec::PetscVec)
  VecAssemblyEnd(vec)
end

function assembly_end(vec::AbstractVector)
end


function getLocalIndices(vec::PetscVec)

  low, high = VecGetOwnershipRange(vec)
  return Int(low):Int(high-1)
end


import Base.maximum
function maximum(x::PetscVec)
  mval, idx = VecMax(x)
  return mval
end

import Base.minimum
function min(x::PetscVec)
  mval, idx = VecMin(x)
  return mval
end

import Base.copy
function copy(vec::PetscVec)
  b2 = VecDuplicate(vec)
  VecCopy(vec, b2)

  return b2
end

#TODO: length_local, length_global, size_local, size_global, fill_zero
#      +, .*, -, (check for in-place versions), in-place transpose, copy!
#      norm, fill!, maximum, minimum


#TODO: extend blas routines

import Base.scale!
function scale!(vec::PetscVec, a::Number)
  _a = PetscScalar(a)
  VecScale(vec, _a)
end

import Base.dot
"""
  Dot product where the *first* vector is conjugated.  This is is the reverse
  of VecDot, where the *second* vector is conjugated
"""
function dot(x::PetscVec, y::PetscVec)
  VecDot(y, x)
end

import Base.sum
function sum(x::PetscVec)
  VecSum(x)
end


