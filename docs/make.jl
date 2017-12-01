using Documenter, PETSc

makedocs(
  format = :html,
  sitename = "Petsc.jl Documentation",
  pages = [
    "index.md",
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

