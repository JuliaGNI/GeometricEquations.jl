
using SafeTestsets

@safetestset "General Equation Functionality                                                  " begin include("general_equations_tests.jl") end
@safetestset "Differential Equations                                                          " begin include("differential_equations_tests.jl") end
@safetestset "Differential-Algebraic Equations                                                " begin include("differential_algebraic_equations.jl") end
@safetestset "Stochastic Equations                                                            " begin include("stochastic_equations_tests.jl") end
