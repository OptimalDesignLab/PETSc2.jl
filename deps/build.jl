
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MPI")
  Pkg.clone("https://github.com/JaredCrean2/MPI.jl.git")
  Pkg.build("MPI")
end

if !haskey(pkg_dict, "ArrayViews")
  Pkg.clone("https://github.com/JuliaLang/ArrayViews.jl.git")
  start_dir = pwd()
  cd(Pkg.dir("ArrayViews"))
  run(`git checkout 715d637fbdd0f6dca6a03d1f2d7c8d248b4abde8`)
  cd (start_dir)
end

if !(haskey(pkg_dict, "FactCheck"))
  Pkg.clone("FactCheck")
  start_dir = pwd()
  cd(Pkg.dir("FactCheck"))
  run(`git checkout 66fc925bd99b50aacfa0f9cd2c324f80ec4947f3`)
  cd(start_dir)
  Pkg.build("FactCheck")
end

if !haskey(ENV, "PETSC_DIR")  && !haskey(ENV, "PETSC_ARCH")
 run(`./install_petsc.sh`)
end 
