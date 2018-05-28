# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.6")



#=
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MPI")
  Pkg.clone("https://github.com/JaredCrean2/MPI.jl.git")
  Pkg.build("MPI")
end

if !haskey(pkg_dict, "ArrayViews")
  Pkg.clone("https://github.com/JuliaLang/ArrayViews.jl.git")
  start_dir = pwd()
  cd(Pkg.dir("ArrayViews"))
  run(`git checkout 93e80390aeedb1dbcd90281b6dff7f760f430bc8`)
  cd(start_dir)
end

if !(haskey(pkg_dict, "FactCheck"))
  Pkg.clone("FactCheck")
  start_dir = pwd()
  cd(Pkg.dir("FactCheck"))
  run(`git checkout e3739d5fdf0e54bc1e74957c060c693cd8ce9cd6`)
  cd(start_dir)
  Pkg.build("FactCheck")
end
=#


global const ARRAYVIEWS_URL = "http://github.com/JuliaArrays/ArrayViews.jl.git"
global const ARRAYVIEWS_VER = "master"

pkg_dict = PkgFix.installed()

if !haskey(pkg_dict, "ArrayViews")
  PkgFix.add(ARRYAVIEWS_URL, branch_ish=ARRAYVIEWS_VER)
else
  PkgFix.checkout("ArrayViews", ARRAYVIEWS_VER)
end




if !haskey(ENV, "PETSC_DIR")  && !haskey(ENV, "PETSC_ARCH")
 run(`./install_petsc.sh`)
end 
