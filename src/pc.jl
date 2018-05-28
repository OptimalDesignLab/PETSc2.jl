export PC, KSPGetPC, PCSetType, PCGetType, PCFactorSetUseInPlace, PCFactorGetUseInPlace, PCSetReusePreconditioner, PCGetReusePreconditioner, PCFactorSetAllowDiagonalFill, PCFactorGetAllowDiagonalFill, PCFactorSetLevels, PCFactorGetLevels, PCSetReusePreconditioner, PCGetReusePreconditioner, PCBJacobiGetSubKSP, PCFactorSetFill, PCJacobiSetType, PCJacobiGetType

#TODO: remove PCShell prefix?
export PCShellSetApply, PCShellSetApplyTranspose, PCShellSetSetUp,
       PCShellSetContext, PCShellGetContext

# developer PC interface
export PCApply, PCSetUp, PCApplyTranspose, PCApplyTransposeExists

export KSPSetPC, KSPGetPC

# preconditioner contex
# the KSP object creates the PC contex, so we don't provide a constructor
"""
  Petsc PC object
"""
type PC
  pobj::Ptr{Void}
end

"""
  Constructor

  **Inputs**

   * comm: MPI communicator

  **Outputs**

   * PC object
"""
function PC(comm::MPI.Comm)
  arg2 = Ref{Ptr{Void}}() 
  ccall((:PCCreate,petsc),PetscErrorCode,(comm_type,Ptr{Ptr{Void}}),comm,arg2)

  return PC(arg2[])
end

"""
  Free a Petsc PC object.  Safe to call multiple times
"""
function PetscDestroy(arg1::PC)
  if arg1.pobj != C_NULL
    ccall((:PCDestroy,petsc),PetscErrorCode,(Ptr{Ptr{Void}},),&arg1.pobj)
  end
end

"""
  PCSetFromOptions
"""
function SetFromOptions(arg1::PC)
    ccall((:PCSetFromOptions,petsc),PetscErrorCode,(Ptr{Void},),arg1.pobj)
end


"""
  KSPGetPC

  **Inputs**

   * ksp: KSP object

  **Output**

   * PC: pc object
"""
function KSPGetPC(ksp::KSP)
    arr = Array(Ptr{Void}, 1)
    ccall((:KSPGetPC,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),ksp.pobj, arr)
    return PC(arr[1])
end

"""
  PCSetType
"""
function PCSetType(pc::PC, pctype::PCType)
    ccall((:PCSetType,petsc),PetscErrorCode,(Ptr{Void},Cstring), pc.pobj, pctype)
end

"""
  KSPSetPC
"""
function KSPSetPC(ksp::KSP, pc::PC)
    ccall((:KSPSetPC,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),ksp.pobj, pc.pobj)
    return nothing
end


"""
  PCGetType

  **Inputs**

   * pc: PC object

  **Outputs**

   * string contianing the PC type
"""
function PCGetType(pc::PC)
    arr = Array(Ptr{UInt8}, 1)
    ccall((:PCGetType,petsc),PetscErrorCode,(Ptr{Void},Ptr{Ptr{UInt8}}), pc.pobj, arr)
    return unsafe_string(arr[1])
end

"""
  PCFactorSetUseInPlace
"""
function PCFactorSetUseInPlace(pc::PC, arg2::PetscBool)
    ccall((:PCFactorSetUseInPlace,petsc),PetscErrorCode,(Ptr{Void},PetscBool),pc.pobj, arg2)
end

"""
  PCFactorGetUseInPlace

  **Inputs**

   * pc: PC object

  **Outputs**

   * PetscBool
"""
function PCFactorGetUseInPlace(pc::PC)
    arr = Array(PetscBool, 1)
    ccall((:PCFactorGetUseInPlace,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end


"""
  PCSetReusePreconditioner
"""
function PCSetReusePreconditioner(pc::PC,arg2::PetscBool)
    ccall((:PCSetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void},PetscBool), pc.pobj, arg2)
end

"""
  PCGetReusePreconditioner

  **Inputs**

   * pc: PC object

  **Outputs**

   * PetscBool
"""
function PCGetReusePreconditioner(pc::PC)
    arr = Array(PetscBool, 1)
    ccall((:PCGetReusePreconditioner,petsc),PetscErrorCode,(Ptr{Void}, Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end

"""
  PCFactorSetAllowDiagonalFill
"""
function PCFactorSetAllowDiagonalFill(pc::PC,arg2::PetscBool)
    ccall((:PCFactorSetAllowDiagonalFill,petsc),PetscErrorCode,(Ptr{Void},PetscBool), pc.pobj, arg2)
end

"""
  PCFactorGetAllowDiagonalFill

  **Inputs**

   * pc: PC object

  **Outputs**

   * PetscBool
"""
function PCFactorGetAllowDiagonalFill(pc::PC)
   arr = Array(PetscBool, 1)
    ccall((:PCFactorGetAllowDiagonalFill,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscBool}), pc.pobj, arr)
    return arr[1]
end

"""
  PCFactorSetLevels
"""
function PCFactorSetLevels(pc::PC,arg2::PetscInt)
    ccall((:PCFactorSetLevels,petsc),PetscErrorCode,(Ptr{Void}, PetscInt), pc.pobj, arg2)
end

"""
  PCFactorGetLevels

  **Inputs**

   * pc: PC object

  **Outputs**

   * PetscInt
"""
function PCFactorGetLevels(pc::PC)
    arr = Array(PetscInt, 1)
    ccall((:PCFactorGetLevels,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt}),pc.pobj, arr)
    return arr[1]
end

"""
  PCBJacobiGetSubKSP

  **Inputs**

   * pc: PC object

  **Outputs**

   * n_local: number of local KSP object
   * first_local: global number of first KSP object on this block
   * ksp_arr: array of KSP object
"""
function PCBJacobiGetSubKSP(pc::PC)
    n_local_arr = Array(PetscInt, 1)
    first_local = Array(PetscInt, 1)
    ksp_ptrarr = Array(Ptr{Ptr{Void}}, 1)
    ccall((:PCBJacobiGetSubKSP,petsc),PetscErrorCode,(Ptr{Void},Ptr{PetscInt},Ptr{PetscInt},Ptr{Ptr{Ptr{Void}}}), pc.pobj, n_local_arr, first_local, ksp_ptrarr)

    n_local = n_local_arr[1]
    ksp_ptrarr2 = unsafe_wrap(Array, ksp_ptrarr[1], n_local)

    ksp_arr = Array(KSP, n_local)
    for i=1:n_local
      ksp_arr[i] = KSP(ksp_ptrarr2[i])
    end

    return n_local, first_local[1], ksp_arr
end

"""
  PCFactorSetFill
"""
function PCFactorSetFill(pc::PC, fill::PetscReal)
    ccall((:PCFactorSetFill,petsc),PetscErrorCode,(Ptr{Void}, PetscReal), pc.pobj, fill)
end

"""
  PCJacobiaSetType
"""
function PCJacobiSetType(pc::PC, jacobitype::PCJacobiType)
    ccall((:PCJacobiSetType,petsc),PetscErrorCode,(Ptr{Void}, PCJacobiType), pc.pobj, jacobitype)
end

"""
  PCJacobiGetType

  **Inputs**

   * pc: PC object

  **Outputs**

   * string containing the type
"""
function PCJacobiGetType(pc::PC)
    arr = Array(PCJacobiType, 1)
    ccall((:PCJacobiGetType,petsc),PetscErrorCode,(Ptr{Void},Ptr{PCJacobiType}),pc.pobj, arr)
    return arr[1]
end

"""
  PCShellSetApply
"""
function PCShellSetApply(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApply,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1.pobj,arg2)
end

"""
  PCShellSetApplyTranspose
"""
function PCShellSetApplyTranspose(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetApplyTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1.pobj,arg2)
end

"""
  PCShellSetSetUp
"""
function PCShellSetSetUp(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetSetUp,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1.pobj,arg2)
end

"""
  PCShellSetContext
"""
function PCShellSetContext(arg1::PC,arg2::Ptr{Void})
    ccall((:PCShellSetContext,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void}),arg1.pobj,arg2)
end

"""
  PCShellGetContext

  **Inputs**

   * pc: PC object

  **Outputs**

   * Ptr{Void}.  Users shoudl call unsafe_pointer_to_objref() on it to get the
                 Julia object back
"""
function PCShellGetContext(arg1::PC)
    arg2 = Ref{Ptr{Void}}()
    ccall((:PCShellGetContext,petsc),PetscErrorCode,(Ptr{Void},Ref{Ptr{Void}}),arg1.pobj,arg2)

    return arg2[]
end

"""
  PCApply
"""
function PCApply(arg1::PC,arg2::PetscVec,arg3::PetscVec)
    ccall((:PCApply,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),arg1.pobj ,arg2.pobj, arg3.pobj)
end

"""
  PCApplyTranspose
"""
function PCApplyTranspose(arg1::PC,arg2::PetscVec,arg3::PetscVec)
    ccall((:PCApplyTranspose,petsc),PetscErrorCode,(Ptr{Void},Ptr{Void},Ptr{Void}),arg1.pobj ,arg2.pobj, arg3.pobj)
end

"""
  PCApplyTransposeExists

  **Inputs**

   * pc: PC object

  **Outputs**

   * Bool
"""
function PCApplyTransposeExists(arg1::PC)
    arg2 = Ref{PetscBool}()
    ccall((:PCApplyTransposeExists,petsc),PetscErrorCode,(PC,Ptr{PetscBool}),arg1,arg2)

    return arg2[] == 1
end

"""
  PCSetUp
"""
function PCSetUp(arg1::PC)
    ccall((:PCSetUp,petsc),PetscErrorCode,(Ptr{Void},),arg1.pobj)
end


