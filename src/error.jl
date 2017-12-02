# have Petsc report errors to Julia
"""
  Error handler registered with PETSc.  This function prints the information
  Petsc supplies about the error, a Julia stack trace, and then invokes the
  `PetscTraceBackErrorHandler` to print the Petsc stack trace.
"""
function error_handler(comm::comm_type, line::Cint, func::Ptr{UInt8}, file::Ptr{UInt8}, n::PetscErrorCode, p::PetscErrorType, mess::Ptr{UInt8}, ctx::Ptr{Void})
# receives the error call from Petsc

  func_string = bytestring(func)
  file_string = bytestring(file)

  if p == PETSC_ERROR_INITIAL
    p_string = "Initial Petsc Error"
  elseif p == PETSC_ERROR_REPEAT
    p_string = "Repeat Petsc Error"
  elseif  p == PETSC_ERROR_IN_CXX
    p_string == "CXX Petsc Error"
  else
    pstring = "Unknown PetscErrorType"
  end

  mess_string = bytestring(mess)
  print_with_color(:red, STDERR, string("\n### ERROR: PETSc Internal Error ###\n"))
  print_with_color(:red, STDERR, string("Error in function ", func_string, ", file ", file_string, "\n"))
  print_with_color(:red, STDERR, string("Error number: ", n, ", Error type: ", p_string, "\n"))
  print_with_color(:red, STDERR, string("Error message: \n", mess_string, "\n"))
  print_with_color(:red, STDERR, string("### Error Message Finished ###\n"))

  bt = backtrace()
  s = sprint(io->Base.show_backtrace(io, bt))
  print_with_color(:red, STDERR, string("Julia backtrace: ", s))
  print(STDERR, "\n\n")

  # invoking Petsc error handler
  PetscTraceBackErrorHandler(comm, line, func, file, n, p, mess, ctx)

  println(STDERR, "finished invoking Petsc error handler")

  return PetscErrorCode(0)

end

function PetscTraceBackErrorHandler(arg1::comm_type,arg2::Cint,arg3::Ptr{UInt8},
           arg4::Ptr{UInt8},arg5::PetscErrorCode,arg6::PetscErrorType,
           arg7::Ptr{UInt8},arg8::Ptr{Void})

    ccall((:PetscTraceBackErrorHandler,petsc),PetscErrorCode,
      (comm_type,Cint,Ptr{UInt8},Ptr{UInt8},PetscErrorCode,PetscErrorType,Ptr{UInt8},Ptr{Void}),
      arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)
end


# tell Petsc abut the Error handler
cfunc = cfunction(error_handler, PetscErrorCode, (comm_type, Int32, Ptr{UInt8}, Ptr{UInt8}, PetscErrorCode, PetscErrorType, Ptr{UInt8}, Ptr{Void}) )
ctx = C_NULL
ierr = ccall( (:PetscPushErrorHandler, petsc), PetscErrorCode, (Ptr{Void}, Ptr{Void}), cfunc, ctx)


function chkerrq(i::PetscErrorCode)
  if i != 0
    error("Petsc Error Code: $i")
  end

  return nothing
end


