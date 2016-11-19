using Documenter, Networks

makedocs(
  format = :html,
  sitename = "Networks",
  pages = [
    "index.md"
  ]
)

deploydocs(
  repo = "github.com/mpichl87/Networks.jl.git",
  target = "build",
  latest = "master",
  julia = "0.5",
  deps = nothing,
  make = nothing
)
