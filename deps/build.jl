
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MPI")
  Pkg.clone("MPI")
  Pkg.build("MPI")
end

if !(haskey(pkg_dict, "FactCheck"))
  Pkg.clone("FactCheck")
  Pkg.checkout("FactCheck", "v0.4.2")
end

if !haskey(ENV, "PETSC_DIR")  && !haskey(ENV, "PETSC_ARCH")
 run(`./install_petsc.sh`)
end 
