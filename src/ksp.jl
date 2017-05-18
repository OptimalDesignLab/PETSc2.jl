
export KSP, SetOperators, SetFromOptions, KSPSolve,  SetUp

export GetConvergedReason, PetscView, SetType, GetType, SetTolerances, GetTolerances, SetInitialGuessNonzero, GetInitialGuessNonzero, GetResidualNorm

type KSP
  pobj::Ptr{Void}

  function KSP(comm::MPI_Comm)  # constructor for KSP owned by the user
      ptr = Array(Ptr{Void}, 1)
      ierr = ccall((:KSPCreate,petsc),PetscErrorCode,(comm_type, Ptr{Void}),comm, ptr)
      @assert(ierr == 0)
      obj = new(ptr[1])
#      finalizer(obj, PetscDestroy)
      return obj
  end

  function KSP(ptr::Ptr{Void})  # constructor for KSP *not* owned by the user
    return new(ptr)
  end
end


function PetscDestroy(ksp::KSP)
    err = ccall((:KSPDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),&ksp.pobj)
    println("KSPDestroy called")
#    sleep(5)
end

function PetscView(ksp::KSP, viewer)
    ccall((:KSPView,petsc),PetscErrorCode,(Ptr{Void}, Int64), ksp.pobj, viewer)
end



function SetOperators(ksp::KSP,Amat::PetscMat,Pmat::PetscMat)
   err = ccall((:KSPSetOperators,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, Amat.pobj, Pmat.pobj)
end


function SetFromOptions(ksp::KSP)
    ccall((:KSPSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},),ksp.pobj)
end

function KSPSolve(ksp::KSP, b::PetscVec, x::PetscVec)
    err = ccall((:KSPSolve,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}), ksp.pobj, b.pobj, x.pobj)
end

function SetUp(ksp::KSP)
    err = ccall((:KSPSetUp,petsc),PetscErrorCode,(Ptr{Void},), ksp.pobj)
end

function GetConvergedReason(ksp::KSP)
    reason = Array(KSPConvergedReason, 1)
    ccall((:KSPGetConvergedReason,petsc),PetscErrorCode,(Ptr{Void},Ptr{KSPConvergedReason}), ksp.pobj, reason)

    return reason[1]
end


function SetType(ksp::KSP, ksptype::KSPType)
    ccall((:KSPSetType,petsc),PetscErrorCode,(Ptr{Void},Cstring), ksp.pobj, ksptype)
end

function GetType(ksp::KSP)
    ksptype = Array(Ptr{UInt8}, 1)
    ccall((:KSPGetType,petsc),PetscErrorCode,(Ptr{Void}, Ptr{Ptr{UInt8}}), ksp.pobj, ksptype)

    return bytestring(ksptype[1])
end


function SetTolerances(ksp::KSP,rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, maxits::PetscInt)

    ccall((:KSPSetTolerances,petsc),PetscErrorCode,(Ptr{Void}, PetscReal, PetscReal, PetscReal, PetscInt), ksp.pobj, rtol, abstol, dtol, maxits)
end

function GetTolerances(ksp::KSP)
    rtol = Array(PetscReal,1)
    abstol = Array(PetscReal,1)
    dtol = Array(PetscReal,1)
    maxits = Array(PetscInt, 1)

    ccall((:KSPGetTolerances,petsc),PetscErrorCode,(KSP,Ptr{PetscReal},Ptr{PetscReal},Ptr{PetscReal},Ptr{PetscInt}), ksp, rtol, abstol, dtol, maxits)

    return rtol[1], abstol[1], dtol[1], maxits[1]
end

function SetInitialGuessNonzero(ksp::KSP, flg::PetscBool)
    ccall((:KSPSetInitialGuessNonzero,petsc),PetscErrorCode,(Ptr{Void},PetscBool), ksp.pobj, flg)
end

function GetInitialGuessNonzero(ksp::KSP)
    flg_arr = Array(PetscBool, 1)
    ccall((:KSPGetInitialGuessNonzero,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscBool}), ksp.pobj, flg_arr)
    return flg_arr[1]
end


function GetResidualNorm(ksp::KSP)
   rnorm = Array(PetscReal, 1)
   ccall((:KSPGetResidualNorm,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscReal}),ksp.pobj, rnorm)

    return rnorm[1]
end





### new function
# not tested
function SetReusePreconditioner(arg1::KSP,arg2::PetscBool)
    ccall((:KSPSetReusePreconditioner,petsc),PetscErrorCode,(KSP,PetscBool),arg1,arg2)
end

