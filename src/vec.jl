export PetscVec, VecSetType, VecSetValues, VecAssemblyBegin, VecAssemblyEnd, VecSetSizes, VecGetSize, VecGetLocalSize, VecNorm, VecGetValues, VecGetOwnershipRange, VecGetArray, VecRestoreArray, VecGetArrayRead, VecRestoreArrayRead, VecSet, VecSqrtAbs, VecLog, VecExp, VecAbs, VecMax, VecMin, VecCopy, VecDuplicate, VecAXPY, VecAXPBY, VecAYPX, VecWAXPY, VecMAXPY, VecAXPBYPCZ, VecScale, VecDot, VecTDot, VecSum, VecSwap, VecReciprocal, VecShift, VecPointwiseMult, VecPointwiseDivide

export getLocalIndices

"""
Petsc Vector type.

Not a subtype of `AbstractArray` because Petsc vectors do not conform to that
API.
"""
type PetscVec
  pobj::Ptr{Void}

  """
    Constructor

    **Inputs**

     * comm: MPI Communicator
  """
  function PetscVec(comm::MPI_Comm)
#    comm = PETSC_COMM_SELF();
    vec = Array{Ptr{Void}}(1)
    err = ccall(( :VecCreate, libpetsclocation ), PetscErrorCode,(comm_type,Ptr{Void}),comm,vec);
    vec = new(vec[1])
#    finalizer(vec,PetscDestroy)
    # does not seem to be called immediately when vec is no longer visible, is it called later during garbage collection? - yes
    return vec

  end

  """
    Constructor from existing pointer to a Petsc Vec

    **Inputs**

     pobj: a Ptr{Void} to an existing Petsc object
  """
  function PetscVec(pobj::Ptr{Void})  # default constructor
    return new(pobj)
  end
 
end

"""
Union{AbstractVector, PetscVec}
"""
const AllVectors = Union{AbstractVector, PetscVec}

"""
  Free a Petsc vec.  Safe to call multiple times
"""
function PetscDestroy(vec::PetscVec)
  if (vec.pobj != 0)
    err = ccall(( :VecDestroy, libpetsclocation), PetscErrorCode, (Ptr{Ptr{Void}},), &vec.pobj);
  end
  vec.pobj = 0  # unnecessary? vec no longer has any references to it
#    println("VecDestroy called")
end

"""
  VecSetType
"""
function VecSetType(vec::PetscVec,name)
  err = ccall((:VecSetType,  libpetsclocation), PetscErrorCode, (Ptr{Void}, Cstring), vec.pobj,name);
end

  #TODO: make this better: use VecFromArray
function PetscVec(array::Array{PetscScalar})
  vec = PetscVec()
  err = ccall(( :VecSetType,  libpetsclocation), PetscErrorCode,(Ptr{Void},Cstring), vec.pobj,"seq");
  err = ccall( (:VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,length(array),length(array));
  # want a PetscInt array so build it ourselves
  idx = Array{PetscInt}(length(array));
  for i=1:length(array);  idx[i] = i-1;  end
  err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode,(Ptr{Void},PetscInt, Ptr{PetscInt},Ptr{PetscScalar},Int32), vec.pobj,length(idx),idx,array,INSERT_VALUES);
  err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode,(Ptr{Void},), vec.pobj);
  err = ccall( ( :VecAssemblyEnd,  libpetsclocation), PetscErrorCode, (Ptr{Void},), vec.pobj);
  return vec
end

"""
   VecSetValues
"""
function VecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar},flag::Integer)

  err = ccall( ( :VecSetValues,  libpetsclocation), PetscErrorCode, 
               (Ptr{Void},PetscInt ,Ptr{PetscInt},Ptr{PetscScalar},Int32), 
               vec.pobj,length(idx), idx,array,flag);
  return err
end

"""
  VecSetValues method that implicitly uses `INSERT_VALUES`
"""
function VecSetValues(vec::PetscVec,idx::Array{PetscInt},array::Array{PetscScalar})
  VecSetValues(vec,idx,array,INSERT_VALUES)
end

#=
function VecSetValues(vec::PetscVec,array::Array{PetscScalar})
  idx = Array{PetscInt}(ength(array))
  for i=1:length(array);  idx[i] = i-1;  end
  VecSetValues(vec,idx,array,INSERT_VALUES)
end
=#

"""
  VecAssemblyBegin
"""
function VecAssemblyBegin(obj::PetscVec)
  err = ccall( ( :VecAssemblyBegin,  libpetsclocation), PetscErrorCode, (Ptr{Void},), obj.pobj);
end

"""
  VecAssemblyEnd
"""
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

"""
  VecSetSizes
"""
function VecSetSizes(vec::PetscVec,n::Integer, N::Integer)
  err = ccall( ( :VecSetSizes,  libpetsclocation), PetscErrorCode, (Ptr{Void}, PetscInt, PetscInt), vec.pobj,n,N);
end

"""
  PetscView for Petsc vector
"""
function PetscView(obj::PetscVec,viewer)
  err = ccall( ( :VecView,  libpetsclocation), PetscErrorCode, (Ptr{Void},Int64),obj.pobj,0);
end

"""
  VeGetSize

  **Inputs**

   * vec: a Petsc vector

  **Outputs**

   * the size
"""
function VecGetSize(obj::PetscVec)
  n = Array{PetscInt}(1)
  err = ccall( ( :VecGetSize,  libpetsclocation), PetscErrorCode, (Ptr{Void},Ptr{PetscInt}), obj.pobj,n);
  return n[1]
end

"""
  VecGetLocalSize

  **Inputs**

   * vec: a Petsc vector

  **Outputs**

   * the size

"""
function VecGetLocalSize(arg1::PetscVec)
  arg2 = Ref{PetscInt}()
  ccall((:VecGetLocalSize,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt}),arg1.pobj,arg2)
  return arg2[]
end



"""
  VecNorm

  **Inputs**

   * obj: Petsc vector
   * normtype: the Petsc enum for the norm type

  **Output**

   * the norm value
"""
function VecNorm(obj::PetscVec,normtype::Integer)
  n = Array{PetscReal}()
  err = ccall( ( :VecNorm,  libpetsclocation), PetscScalar, (Ptr{Void},Int32,Ptr{PetscReal}), obj.pobj,normtype, n);
  return n[1]
end

#=
function VecNorm(obj::PetscVec)
  return VecNorm(obj,PETSC_NORM_2)
end
=#
"""
  VecGetValues
"""
function VecGetValues(vec::PetscVec, ni::Integer, ix::AbstractArray{PetscInt,1}, y::AbstractArray{PetscScalar,1})

     # need indices to be PetscInt
#     ix_local = Array{PetscInt}(ni)
#     for i=1:ni
#       ix_local[i] = ix[i]
#     end

    err = ccall((:VecGetValues,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscInt},Ptr{PetscScalar}),vec.pobj, ni, ix, y)

    return nothing
end

"""
  VecGetValues with length of idx inferred from idx

  **Inputs**

   * vec: the Petsc vector
   * idx: array of PetscInt indices
   * y: array of PetscScalar values
"""
function VecGetValues(vec::PetscVec,idx::AbstractArray{PetscInt,1}, y::AbstractArray{PetscScalar, 1})

  VecGetValues(vec, length(idx), idx, y)
end

"""
  VecGetOwnershipRange

  **Inputs**

   * vec: Petsc vector

  **Outputs**

   * low: lowest index (zero-based) that is owned
   * high: highest index + 1 that is owned
"""
function VecGetOwnershipRange(vec::PetscVec)
    low = Array{PetscInt}(1)
    high = Array{PetscInt}(1)

    ccall((:VecGetOwnershipRange,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt}),vec.pobj, low, high)

  return low[1], high[1]
end


# new functions
"""
  VecGetArray.  Users must call `VecRestoreArray` when finished.

  **Inputs**

   * vec: the Petsc vector

  **Outputs**

   * arr: a Julia Array{PetscScalar, 1}
"""
function VecGetArray(vec::PetscVec)
# gets a pointer to the data underlying a Petsc vec, turns it into a Julia
# array
# ptr_arr must be passed into PetscVecRestoreArray
    ptr_arr = Array{Ptr{PetscScalar}}(1)
    ccall((:VecGetArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = VecGetOwnershipRange(vec)
    len = last - first
    arr = unsafe_wrap(Array, ptr_arr[1], len)
    return arr
end

"""
  VecRestoreArray.  Users must not access the array after calling this function.

  **Inputs**

   * vec: the PetscVector passed into `VecGetArray`
   * arr: the array returned by `VecGetArray
"""
function VecRestoreArray(vec::PetscVec, arr::Array{PetscScalar, 1})
    ptr_arr = Ref{Ptr{PetscScalar}}(pointer(arr))
    ccall((:VecRestoreArray,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
end


"""
  Similar to `VecGetArray`, but produces the array must not be written to.
"""
function VecGetArrayRead(vec::PetscVec)
    ptr_arr = Array{Ptr{PetscScalar}}(1)
    ccall((:VecGetArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)

    first, last = VecGetOwnershipRange(vec)
    len = last - first
    arr = unsafe_wrap(Array, ptr_arr[1], len)
    return arr
end

"""
  Similar to `VecRestoreArray`, but corresponds to `VecGetArrayRead`
"""
function VecRestoreArrayRead(vec::PetscVec, arr::Array{PetscScalar, 1})

    ptr_arr = Ref{Ptr{PetscScalar}}(pointer(arr))
    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{PetscScalar}}),vec.pobj, ptr_arr)
#    ccall((:VecRestoreArrayRead,petsc),PetscErrorCode,(Vec,Ptr{Ptr{PetscScalar}}),arg1,arg2)
end


"""
  VecSet
"""
function VecSet(vec::PetscVec, val::PetscScalar)
    ccall((:VecSet,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, val)
end

"""
  VecSqrtAbs
"""
function VecSqrtAbs(vec::PetscVec)
    ccall((:VecSqrtAbs,petsc),PetscErrorCode,(Ptr{Void},), vec.pobj)
end

"""
  VecLog
"""
function VecLog(vec::PetscVec)
    ccall((:VecLog,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

"""
  VecExp
"""
function VecExp(vec::PetscVec)
    ccall((:VecExp,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

function VecAbs(vec::PetscVec)
    ccall((:VecAbs,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

"""
  VecMax

  **Inputs**

   * vec: PetscVec

  **Output**

   * r: the maximum value
   * idx: the (zero-based) index of the maximum value
"""
function VecMax(vec::PetscVec)
    r = Array{PetscReal}(1) # max value
    idx = Array{PetscInt}(1)  # index of max value
    ccall((:VecMax,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end

"""
  VecMin. Same interface as `VecMax`
"""
function VecMin(vec::PetscVec)
    r = Array{PetscReal}(1) # min value
    idx = Array{PetscInt}(1)  # index of min value
 
    ccall((:VecMin,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscReal}), vec.pobj, idx, r)

    return r[1], idx[1]
end

"""
  VecReciprocal
"""
function VecReciprocal(vec::PetscVec)
    ccall((:VecReciprocal,petsc),PetscErrorCode,(Ptr{Void},),vec.pobj)
end

"""
  VecShift
"""
function VecShift(vec::PetscVec, a::PetscScalar)
    ccall((:VecShift,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end

"""
  VecPointwiseMult
"""
function VecPointwiseMult(w::PetscVec, x::PetscVec,y::PetscVec)
    ccall((:VecPointwiseMult,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end

"""
  VecPointwiseDivide
"""
function VecPointwiseDivide(w::PetscVec, x::PetscVec, y::PetscVec)
    ccall((:VecPointwiseDivide,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}, Ptr{Void}), w.pobj, x.pobj, y.pobj)
end





"""
  VecCopy
"""
function VecCopy(vec::PetscVec , vec2::PetscVec)
    ccall((:VecCopy,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), vec.pobj, vec2.pobj)
end

"""
  VecDulicate
"""
function VecDuplicate( vec::PetscVec)
    ptr_arr = Array{Ptr{Void}}(1)

    ccall((:VecDuplicate,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Ptr{Void}}), vec.pobj, ptr_arr)

    return PetscVec(ptr_arr[1])
end


# Some vector linear algebra
"""
  VecAXPY
"""
function VecAXPY( vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}), vec1.pobj, a, vec2.pobj)
end

"""
  VecAXPBY
"""
function VecAXPBY(vec1::PetscVec, a::PetscScalar, b::PetscScalar, vec2::PetscVec)
    ccall((:VecAXPBY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,Ptr{Void}),vec1.pobj, a, b, vec2.pobj)
end

"""
  VecMAXPY
"""
function VecMAXPY(vec1::PetscVec, n::Integer, a::AbstractArray{PetscScalar, 1}, x::AbstractArray{Ptr{Void}, 1})
# the vector x must contains the pointers from the PetscVec objects, not the PetscVec objects themselves

    ccall((:VecMAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscInt,Ptr{PetscScalar},Ptr{Ptr{Void}}),vec1.pobj, n, a, x)
end

"""
  VecAYPX
"""
function VecAYPX(vec1::PetscVec, a::PetscScalar, vec2::PetscVec)
    ccall((:VecAYPX,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void}),vec1.pobj, a, vec2.pobj)
end

"""
  VecWAXPY
"""
function VecWAXPY(w::PetscVec, a::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecWAXPY,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,Ptr{Void},Ptr{Void}), w.pobj, a, x.pobj, y.pobj)
end

"""
  VecAXPBYPCZ
"""
function VecAXPBYPCZ(z::PetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::PetscVec, y::PetscVec)
    ccall((:VecAXPBYPCZ,petsc),PetscErrorCode,(Ptr{Void},PetscScalar,PetscScalar,PetscScalar,Ptr{Void},Ptr{Void}), z.pobj, alpha,beta, gamma, x.pobj, y.pobj)
end

"""
  VecScale
"""
function VecScale(vec::PetscVec, a::PetscScalar)
    ccall((:VecScale,petsc),PetscErrorCode,(Ptr{Void},PetscScalar), vec.pobj, a)
end

"""
  VecDot
"""
function VecDot(x::PetscVec, y::PetscVec)
    r = Array{PetscScalar}(1)
    err = ccall((:VecDot,petsc),PetscErrorCode,( Ptr{Void}, Ptr{Void}, Ptr{PetscScalar}), x.pobj, y.pobj, r)

    println("PetscVecDot error code = ", err)
    return r[1]
end

"""
  VecTDot
"""
function VecTDot(x::PetscVec, y::PetscVec)
    r = Array{PetscScalar}(1)
    ccall((:VecTDot,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void},Ptr{PetscScalar}), x.pobj, y.pobj, r)

    return r[1]
end

"""
  VecSum
"""
function VecSum(vec::PetscVec)
    r = Array{PetscScalar}(1)
    ccall((:VecSum,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscScalar}),vec.pobj, r)
    return r[1]
end

"""
  VecSwap
"""
function VecSwap(x::PetscVec, y::PetscVec)
    ccall((:VecSwap,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Void}), x.pobj, y.pobj)

end

include("vec_interface.jl")
