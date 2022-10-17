using GeometricEquations
using Documenter

DocMeta.setdocmeta!(GeometricEquations, :DocTestSetup, :(using GeometricEquations); recursive=true)

makedocs(;
    modules=[GeometricEquations],
    authors="Michael Kraus",
    repo="https://github.com/JuliaGNI/GeometricEquations.jl/blob/{commit}{path}#{line}",
    sitename="GeometricEquations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGNI.github.io/GeometricEquations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Equations" => "equations.md",
        "Problems" => "problems.md",
        "Tests" => "tests.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/GeometricEquations.jl",
)
