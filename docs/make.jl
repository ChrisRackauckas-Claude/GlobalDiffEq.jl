using Documenter, GlobalDiffEq

DocMeta.setdocmeta!(
    GlobalDiffEq,
    :DocTestSetup,
    :(using GlobalDiffEq, OrdinaryDiffEq, SciMLBase);
    recursive = true,
)

makedocs(;
    sitename = "GlobalDiffEq.jl",
    modules = [GlobalDiffEq],
    checkdocs = :exports,
    doctest = true,
    pages = [
        "Home" => "index.md",
        "Public API" => "api.md",
    ],
)

deploydocs(repo = "github.com/SciML/GlobalDiffEq.jl.git", devbranch = "master")
