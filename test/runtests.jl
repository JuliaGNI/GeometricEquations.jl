
using SafeTestsets

@safetestset "Utility Functions                                                               " begin include("utils_tests.jl") end
@safetestset "Abstract Equations                                                              " begin include("geometric_equations_tests.jl") end
@safetestset "Ordinary Differential Equations                                                 " begin include("ordinary_differential_equations_tests.jl") end
@safetestset "Differential Algebraic Equations                                                " begin include("differential_algebraic_equations_tests.jl") end
@safetestset "Stochastic Differential Equations                                               " begin include("stochastic_differential_equations_tests.jl") end
