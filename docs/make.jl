using Documenter, Networks

makedocs(
  format = :html,
  sitename = "Networks",
  pages = [
    "index.md"
  ]
)

deploydocs(
  repo = "github.com/mpichl87/Networks.jl.git"
)
