# Interacting with the Petsc Options Database
export PetscOptionsSetValue, PetscOptionsClearValue, PetscOptionsView, PetscSetOptions, PetscClearOptions


typealias PetscOptions Ptr{Void}

function PetscOptionsSetValue(arg1::AbstractString,arg2::AbstractString, arg3=C_NULL)
    ccall((:PetscOptionsSetValue,petsc),PetscErrorCode,(PetscOptions, Cstring, Cstring),arg3, arg1,arg2)
end

function PetscOptionsView(arg1::PetscViewer=C_NULL, arg2=C_NULL)
    ccall((:PetscOptionsView,petsc),PetscErrorCode,(PetscOptions, PetscViewer,),arg2, arg1)
end

function PetscOptionsClearValue(arg1::AbstractString, arg3=C_NULL)
    ccall((:PetscOptionsClearValue,petsc),PetscErrorCode,(PetscOptions, Cstring,),arg3, arg1)
end


"""
  Convenience wrapper for using a dictionary to set options
"""
function PetscSetOptions(opts::Dict)

  for (key, value) in opts
    PetscOptionsSetValue(key, value)
  end

end

function PetscClearOptions(opts::Dict)

  for key in keys(opts)
    PetscOptionsClearValue(key)
  end

end


