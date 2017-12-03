
export KSP, SetOperators, SetFromOptions, KSPSolve,  SetUp, KSPSolveTranspose, KSP_NULL


export GetConvergedReason, PetscView, SetType, GetType, SetTolerances, GetTolerances, SetInitialGuessNonzero, GetInitialGuessNonzero, GetResidualNorm

"""
  KSP object
"""
type KSP
  pobj::Ptr{Void}

  """
    Constructor

    **Inputs**

     * comm: MPI_Comm

    **Outputs**

     * ksp object
  """
  function KSP(comm::MPI_Comm)  # constructor for KSP owned by the user
      ptr = Array(Ptr{Void}, 1)
      ierr = ccall((:KSPCreate,petsc),PetscErrorCode,(comm_type, Ptr{Void}),comm, ptr)
      @assert(ierr == 0)
      obj = new(ptr[1])
#      finalizer(obj, PetscDestroy)
      return obj
  end

  """
    KSP constructor from pointer to existing KSP object
  """
  function KSP(ptr::Ptr{Void})  # constructor for KSP *not* owned by the user
    return new(ptr)
  end
end

"""
  Null pointer KSP object
"""
global const KSP_NULL = KSP(C_NULL)

"""
  PetscDestroy for KSP object.  Safe to call multiple times.
"""
function PetscDestroy(ksp::KSP)

  if ksp.pobj != C_NULL
    err = ccall((:KSPDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),&ksp.pobj)
  end

  ksp.pobj = C_NULL
end

"""
  PetscView

  **Inputs**

   * ksp: KSP object
   * viewer: PetscViewer, defaults to Petsc stdout
"""
function PetscView(ksp::KSP, viewer::PetscViewer=C_NULL)
    ccall((:KSPView,petsc),PetscErrorCode,(Ptr{Void}, Int64), ksp.pobj, viewer)
end


"""
  KSPSetOperators
"""
function SetOperators(ksp::KSP,Amat::PetscMat,Pmat::PetscMat)
   err = ccall((:KSPSetOperators,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, Amat.pobj, Pmat.pobj)
end

"""
  SetFromOptions
"""
function SetFromOptions(ksp::KSP)
    ccall((:KSPSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},),ksp.pobj)
end

"""
  KSPSolve
"""
function KSPSolve(ksp::KSP, b::PetscVec, x::PetscVec)
    err = ccall((:KSPSolve,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, b.pobj, x.pobj)
end

"""
  KSPSolveTranspose
"""
function KSPSolveTranspose(ksp::KSP, b::PetscVec, x::PetscVec)
    err = ccall((:KSPSolveTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, b.pobj, x.pobj)
end


"""
  KSPSetUp
"""
function SetUp(ksp::KSP)
    err = ccall((:KSPSetUp,petsc),PetscErrorCode,(Ptr{Void},), ksp.pobj)
end

"""
  KSPGetConvergedReason

  **Inputs**

   * ksp: KSP object

  **Outputs**

   * string containing reason
"""
function GetConvergedReason(ksp::KSP)
    reason = Array(KSPConvergedReason, 1)
    ccall((:KSPGetConvergedReason,petsc),PetscErrorCode,(Ptr{Void},Ptr{KSPConvergedReason}), ksp.pobj, reason)

    return reason[1]
end

"""
  KSPSetType
"""
function SetType(ksp::KSP, ksptype::KSPType)
    ccall((:KSPSetType,petsc),PetscErrorCode,(Ptr{Void},Cstring), ksp.pobj, ksptype)
end

"""
  KSPGetType

  **Inputs**

   * KSP

  **Outputs**

   * string containing KSP type
"""
function GetType(ksp::KSP)
    ksptype = Array(Ptr{UInt8}, 1)
    ccall((:KSPGetType,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{UInt8}}), ksp.pobj, ksptype)

    return bytestring(ksptype[1])
end

"""
  KSPSetTolerances
"""
function SetTolerances(ksp::KSP,rtol::Number, abstol::Number, dtol::Number, maxits::Integer)

    ccall((:KSPSetTolerances,petsc),PetscErrorCode,(Ptr{Void}, PetscReal, PetscReal, PetscReal, PetscInt), ksp.pobj, rtol, abstol, dtol, maxits)
end

"""
  KSPGetTolerances

  **Inputs**

   * KSP

  **Outputs**

   * rtol
   * abstol
   * dtol
   * maxits
"""
function GetTolerances(ksp::KSP)
    rtol = Array(PetscReal,1)
    abstol = Array(PetscReal,1)
    dtol = Array(PetscReal,1)
    maxits = Array(PetscInt, 1)

    ccall((:KSPGetTolerances,petsc),PetscErrorCode,(KSP,Ptr{PetscReal},Ptr{PetscReal},Ptr{PetscReal},Ptr{PetscInt}), ksp, rtol, abstol, dtol, maxits)

    return rtol[1], abstol[1], dtol[1], maxits[1]
end

"""
  KSPSetInitialGuessNonzero
"""
function SetInitialGuessNonzero(ksp::KSP, flg::PetscBool)
    ccall((:KSPSetInitialGuessNonzero,petsc),PetscErrorCode,(Ptr{Void},PetscBool), ksp.pobj, flg)
end

"""
  KSPGetInitialGuessNonzero

  **Inputs**

   * KSP

  **Outputs**

   * PetscBool
"""
function GetInitialGuessNonzero(ksp::KSP)
    flg_arr = Array(PetscBool, 1)
    ccall((:KSPGetInitialGuessNonzero,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscBool}), ksp.pobj, flg_arr)
    return flg_arr[1]
end

"""
  KSPGetResidualNorm

  **Inputs**

   * KSP

  **Outputs**

   * PetscReal
"""
function GetResidualNorm(ksp::KSP)
   rnorm = Array(PetscReal, 1)
   ccall((:KSPGetResidualNorm,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscReal}),ksp.pobj, rnorm)

    return rnorm[1]
end





### new function
# not tested

"""
  KSPSetReusePreconditioner
"""
function SetReusePreconditioner(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetReusePreconditioner,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end


