using Documenter, PETSc2

makedocs(
  format = :html,
  sitename = "Petsc2.jl",
  pages = [
    "index.md",
    "build.md",
    "init.md",
    "vec.md",
    "vec_interface.md",
    "mat.md",
    "mat_interface.md",
    "constants.md",
    "ksp.md",
    "pc.md",
    "options.md",
    "error.md"
  ]
)

