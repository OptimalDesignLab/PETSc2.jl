export PetscVec, VecSetType, VecSetValues, VecAssemblyBegin, VecAssemblyEnd, VecSetSizes, VecGetSize, VecNorm, VecGetValues, VecGetOwnershipRange, VecGetArray, VecRestoreArray, VecGetArrayRead, VecRestoreArrayRead, VecSet, VecSqrtAbs, VecLog, VecExp, VecAbs, VecMax, VecMin, VecCopy, VecDuplicate, VecAXPY, VecAXPBY, VecAYPX, VecWAXPY, VecMAXPY, VecAXPBYPCZ, VecScale, VecDot, VecTDot, VecSum, VecSwap, VecReciprocal, VecShift, VecPointwiseMult, VecPointwiseDivide, set_values1!, get_values1!

export getLocalIndices

type PetscVec
  pobj::Ptr{Void}
  function PetscVec(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array(Ptr{Void},1)
    err = ccall(( :VecCreate, libpetsclocation ), PetscErrorCode,(comm_type,Ptr{Void}),comm,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection? - yes
    return vec

  end

  function PetscVec(pobj::Ptr{Void})  # default constructor
    return new(pobj)
  end
 
end

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


  function PetscDestroy(vec::PetscVec)
    if (vec.pobj != 0)
      err = ccall(( :VecDestroy, libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
    end
    vec.pobj = 0  # unnecessary? vec no longer has any references to it
#    println("VecDestroy called")
  end

  function VecSetType(vec::PetscVec,name)
    err = ccall((:VecSetType,  libpetsclocation), PetscErrorCode, (Ptr{Void}, Cstring), vec.pobj,name);
  end

  #TODO: make this better: use VecFromArray
  function PetscVec(array::Array{PetscScalar})
    vec = PetscVec()
    err = ccall(( :VecSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void},Cstring), vec.pobj,"seq");
    err = ccall( (:VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,length(array),length(array));
    # want a PetscInt array so build it ourselves
    idx = Array(PetscInt, length(array));
    for i=1:length(array);  idx[i] = i-1;  end
    err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode,(Ptr{Void},PetscInt, Ptr{PetscInt},Ptr{PetscScalar},Int32), vec.pobj,length(idx),idx,array,PETSC_INSERT_VALUES);
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},), vec.pobj);
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);
    return vec
  end

  function VecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar},flag::Integer)

    err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode, 
                 (Ptr{Void},PetscInt ,Ptr{PetscInt},Ptr{PetscScalar},Int32), 
                 vec.pobj,length(idx), idx,array,flag);
    return err
  end

  function VecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar})
    VecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end

  function VecSetValues(vec::PetscVec,array::Array{PetscScalar})
    idx = Array(PetscInt,length(array))
    for i=1:length(array);  idx[i] = i-1;  end
    VecSetValues(vec,idx,array,PETSC_INSERT_VALUES)
  end


  # 1-based indexing that unifies regular matrices and Petsc matrices
  #----------------------------------------------------------------------------
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


  function VecAssemblyBegin(obj::PetscVec)
    err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode, (Ptr{Void},), obj.pobj);
  end

  function VecAssemblyEnd(obj::PetscVec)
    err = ccall( ( :VecAssemblyEnd,  libpetsclocation), PetscErrorCode,(Ptr{Void},), obj.pobj);
  end

  """
    Convenience function for calling VecAssemblyBegin, and VecAssemblyEnd, in
    one go.
  """
  function VecAssemble(obj::PetscVec)
    VecAssemblyBegin(obj)
    VecAssemblyEnd(obj)
  end

  function VecSetSizes(vec::PetscVec,n::Integer, N::Integer)
    err = ccall( ( :VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,n,N);
  end

  function PetscView(obj::PetscVec,viewer)
    err = ccall( ( :VecView,  libpetsclocation), PetscErrorCode, (Ptr{Void},Int64),obj.pobj,0);
  end

  function VecGetSize(obj::PetscVec)
    n = Array(PetscInt, 1)
    err = ccall( ( :VecGetSize,  libpetsclocation), PetscErrorCode, (Ptr{Void},Ptr{PetscInt}), obj.pobj,n);
    return n[1]
  end


  #TODO: VecGetLocalSize

  function VecNorm(obj::PetscVec,normtype::Integer)
    n = Array(PetscReal,1)
    err = ccall( ( :VecNorm,  libpetsclocation), PetscScalar, (Ptr{Void},Int32,Ptr{PetscReal}), obj.pobj,normtype, n);
    return n[1]
  end

  function VecNorm(obj::PetscVec)
    return VecNorm(obj,PETSC_NORM_2)
  end

function VecGetValues(vec::PetscVec, ni::Integer, ix::AbstractArray{PetscInt,1}, y::AbstractArray{PetscScalar,1})

     # need indices to be PetscInt
#     ix_local = Array(PetscInt, ni)
#     for i=1:ni
#       ix_local[i] = ix[i]
#     end

    err = ccall((:VecGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),vec.pobj, ni, ix, y)

    return nothing
end

function VecGetValues(vec::PetscVec,idx::AbstractArray{PetscInt,1}, y::AbstractArray{PetscScalar, 1})

  VecGetValues(vec, length(idx), idx, y)
end

function VecGetOwnershipRange(vec::PetscVec)
    low = Array(PetscInt, 1)
    high = Array(PetscInt, 1)

    ccall((:VecGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),vec.pobj, low, high)

  return low[1], high[1]
end

function getLocalIndices(vec::PetscVec)

  low, high = VecGetOwnershipRange(vec)
  return Int(low):Int(high-1)
end

#TODO; add function that returns the range of indices, not low:end + 1

# new functions

function VecGetArray(vec::PetscVec)
# gets a pointer to the data underlying a Petsc vec, turns it into a Julia
# array
# ptr_arr must be passed into PetscVecRestoreArray
    ptr_arr = Array(Ptr{PetscScalar}, 1)
    ccall((:VecGetArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = VecGetOwnershipRange(vec)
    len = last - first
    arr = pointer_to_array(ptr_arr[1], len)
    return arr, ptr_arr
end


function VecRestoreArray(vec::PetscVec, ptr_arr::Array{Ptr{PetscScalar}, 1})
    ccall((:VecRestoreArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
end


function VecGetArrayRead(vec::PetscVec)
    ptr_arr = Array(Ptr{PetscScalar}, 1)
    ccall((:VecGetArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = VecGetOwnershipRange(vec)
    len = last - first
    arr = pointer_to_array(ptr_arr[1], len)
    return arr, ptr_arr
end

function VecRestoreArrayRead(vec::PetscVec, ptr_arr::Array{Ptr{PetscScalar}, 1})

    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
#    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end



function VecSet(vec::PetscVec, val::PetscScalar)
    ccall((:VecSet,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, val)
end


function VecSqrtAbs(vec::PetscVec)
    ccall((:VecSqrtAbs,petsc),PetscErrorCode,(Ptr{Void},), vec.pobj)
end

function VecLog(vec::PetscVec)
    ccall((:VecLog,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function VecExp(vec::PetscVec)
    ccall((:VecExp,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function VecAbs(vec::PetscVec)
    ccall((:VecAbs,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end


function VecMax(vec::PetscVec)
    r = Array(PetscReal, 1) # max value
    idx = Array(PetscInt, 1)  # index of max value
    ccall((:VecMax,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end

import Base.maximum
function maximum(x::PetscVec)
  mval, idx = VecMax(x)
  return mval
end

function VecMin(vec::PetscVec)
    r = Array(PetscReal, 1) # min value
    idx = Array(PetscInt, 1)  # index of min value
 
    ccall((:VecMin,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end
import Base.minimum
function min(x::PetscVec)
  mval, idx = VecMin(x)
  return mval
end

function VecReciprocal(vec::PetscVec)
    ccall((:VecReciprocal,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function VecShift(vec::PetscVec, a::PetscScalar)
    ccall((:VecShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end


function VecPointwiseMult(w::PetscVec, x::PetscVec,y::PetscVec)
    ccall((:VecPointwiseMult,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end

function VecPointwiseDivide(w::PetscVec, x::PetscVec, y::PetscVec)
    ccall((:VecPointwiseDivide,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end






function VecCopy(vec::PetscVec , vec2::PetscVec)
    ccall((:VecCopy,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), vec.pobj, vec2.pobj)
end

function VecDuplicate( vec::PetscVec)
    ptr_arr = Array(Ptr{Void}, 1)

    ccall((:VecDuplicate,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Ptr{Void}}), vec.pobj, ptr_arr)

    return PetscVec(ptr_arr[1])
end

import Base.copy
function copy(vec::PetscVec)
  b2 = VecDuplicate(vec)
  VecCopy(vec, b2)

  return b2
end


# Some vector linear algebra
function VecAXPY( vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}), vec1.pobj, a, vec2.pobj)
end

function VecAXPBY(vec1::PetscVec, a::PetscScalar, b::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPBY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,Ptr{Void}),vec1.pobj, a, b, vec2.pobj)
end

function VecMAXPY(vec1::PetscVec, n::Integer, a::AbstractArray{PetscScalar, 1}, x::AbstractArray{Ptr{Void}, 1})
# the vector x must contains the pointers from the PetscVec objects, not the PetscVec objects themselves

    ccall((:VecMAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscScalar},Ptr{Ptr{Void}}),vec1.pobj, n, a, x)
end

function VecAYPX(vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}),vec1.pobj, a, vec2.pobj)
end

function VecWAXPY(w::PetscVec, a::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecWAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},Ptr{Void}), w.pobj, a, x.pobj, y.pobj)
end

function VecAXPBYPCZ(z::PetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecAXPBYPCZ,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,PetscScalar,Ptr{Void},Ptr{Void}), z.pobj, alpha,beta, gamma, x.pobj, y.pobj)
end


function VecScale(vec::PetscVec, a::PetscScalar)
    ccall((:VecScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end

import Base.scale!
function scale!(vec::PetscVec, a::Number)
  _a = PetscScalar(a)
  VecScale(vecc, _a)
end

function VecDot(x::PetscVec, y::PetscVec)
    r = Array(PetscScalar, 1)
    err = ccall((:VecDot,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Void}, Ptr{PetscScalar}), x.pobj, y.pobj, r)

    println("PetscVecDot error code = ", err)
    return r[1]
end

import Base.dot
"""
  Dot product where the *first* vector is conjugated.  This is is the reverse
  of VecDot, where the *second* vector is conjugated
"""
function dot(x::PetscVec, y::PetscVec)
  VecDot(y, x)
end

function VecTDot(x::PetscVec, y::PetscVec)
    r = Array(PetscScalar, 1)
    ccall((:VecTDot,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void},Ptr{PetscScalar}), x.pobj, y.pobj, r)

    return r[1]
end


function VecSum(vec::PetscVec)
    r = Array(PetscScalar, 1)
    ccall((:VecSum,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscScalar}),vec.pobj, r)
    return r[1]
end

import Base.sum
function sum(x::PetscVec)
  VecSum(x)
end

function VecSwap(x::PetscVec, y::PetscVec)
    ccall((:VecSwap,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), x.pobj, y.pobj)

end


