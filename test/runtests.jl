using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "ExplicitImports" include("quality/explicit_imports.jl")
    @safetestset "Utility Functions" include("utils.jl")
    @safetestset "Abstract Equation" include("geometric_equation.jl")
    @safetestset "Ordinary Differential Equations" include("odes/odes.jl")
    @safetestset "Differential Algebraic Equations" include("daes/daes.jl")
    @safetestset "Stochastic Differential Equations" include("sdes/sdes.jl")
    @safetestset "Stochastic Processes" include("sdes/processes.jl")
    @safetestset "Discrete Equations" include("discrete/dele.jl")
    @safetestset "Geometric Problem" include("geometric_problem.jl")
    @safetestset "Equation Problem" include("problems/equation_problem.jl")
    @safetestset "Ensemble Problem" include("problems/ensemble_problem.jl")
    @safetestset "Conversion" include("conversion.jl")
    @safetestset "Test Problems" include("tests/Tests.jl")
end
