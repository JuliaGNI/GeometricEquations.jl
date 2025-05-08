
using SafeTestsets

@safetestset "Utility Functions                                                               " begin include("utils_tests.jl") end
@safetestset "Abstract Equation                                                               " begin include("geometric_equations_tests.jl") end
@safetestset "Ordinary Differential Equations                                                 " begin include("ordinary_differential_equations_tests.jl") end
@safetestset "Differential Algebraic Equations                                                " begin include("differential_algebraic_equations_tests.jl") end
@safetestset "Stochastic Differential Equations                                               " begin include("stochastic_differential_equations_tests.jl") end
@safetestset "Discrete Equations                                                              " begin include("discrete_equations_tests.jl") end
@safetestset "Geometric Problem                                                               " begin include("geometric_problem_tests.jl") end
@safetestset "Equation Problem                                                                " begin include("equation_problem_tests.jl") end
@safetestset "Ensemble Problem                                                                " begin include("ensemble_problem_tests.jl") end
@safetestset "Test Problems                                                                   " begin include("tests_tests.jl") end
