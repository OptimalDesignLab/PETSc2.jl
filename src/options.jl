# Interacting with the Petsc Options Database
export PetscOptionsSetValue, PetscOptionsClearValue, PetscOptionsView, PetscSetOptions, PetscClearOptions

"""
  Typedef of PetscOptions
"""
typealias PetscOptions Ptr{Void}

"""
  PetscOptionsSetValue

  **Inputs**

   * arg1: the key (string)
   * arg2: the value (string)
   * arg3: the PetscOptions object, defaults to the global options database
"""
function PetscOptionsSetValue(arg1::AbstractString,arg2::AbstractString, arg3::PetscOptions=C_NULL)
    ccall((:PetscOptionsSetValue,petsc),PetscErrorCode,(PetscOptions, Cstring, Cstring),arg3, arg1,arg2)
end

"""
  PetscOptionsView

  **Inputs**

   * arg1: a PetscViewer, defaults to Petsc's stdout
   * arg2: the PetscOptions object, defaults to the global options databse
"""
function PetscOptionsView(arg1::PetscViewer=C_NULL, arg2::PetscOptions=C_NULL)
    ccall((:PetscOptionsView,petsc),PetscErrorCode,(PetscOptions, PetscViewer,),arg2, arg1)
end

"""
  PetscOptionsClearValue

  **Inputs**

   * arg1: the key (string)
   * arg2: the PetscOptions object, defaults to the global options databse
"""
function PetscOptionsClearValue(arg1::AbstractString, arg3::PetscOptions=C_NULL)
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

"""
  Convenience wrapper for using a dictionary to clear options (only the keys
  are used).
"""
function PetscClearOptions(opts::Dict)

  for key in keys(opts)
    PetscOptionsClearValue(key)
  end

end


