using GeometricEquations
using Documenter

DocMeta.setdocmeta!(GeometricEquations, :DocTestSetup, :(using GeometricEquations); recursive = true)

makedocs(;
    modules = [GeometricEquations],
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, :parse_error, :setup_block),
    authors = "Michael Kraus",
    sitename = "GeometricEquations.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://JuliaGNI.github.io/GeometricEquations.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Equations" => "equations.md",
        "Problems" => "problems.md",
        "Ensembles" => "ensembles.md",
        "Index" => "reference.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/GeometricEquations.jl",
    devurl = "latest",
    devbranch = "main",
)
