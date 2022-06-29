
using GeometricEquations
using GeometricEquations: arrtype, datatype
using GeometricEquations: check_initial_conditions, check_parameters
using GeometricEquations: _functions, _solutions, _invariants
using Test

@test AbstractEquationODE  <: GeometricEquation
@test AbstractEquationDAE  <: GeometricEquation
@test AbstractEquationSDE  <: GeometricEquation
@test AbstractEquationPODE <: GeometricEquation
@test AbstractEquationPDAE <: GeometricEquation
@test AbstractEquationPSDE <: GeometricEquation


struct TestEquation <: GeometricEquation{Nothing,Nothing,Nothing} end

testeq = TestEquation()

@test functions(testeq) == NamedTuple()
@test solutions(testeq) == NamedTuple()

@test_throws ErrorException functions(testeq, NullParameters())
@test_throws ErrorException solutions(testeq, NullParameters())

@test_throws ErrorException _functions(testeq)
@test_throws ErrorException _solutions(testeq)
@test_throws ErrorException _invariants(testeq)

@test_throws ErrorException _functions(testeq, NullParameters())
@test_throws ErrorException _solutions(testeq, NullParameters())
@test_throws ErrorException _invariants(testeq, NullParameters())

@test_throws ErrorException parameters(testeq)
@test_throws ErrorException periodicity(testeq)

@test_throws ErrorException check_initial_conditions(testeq, NamedTuple())

@test_throws ErrorException check_parameters(testeq, NullParameters())
@test_throws ErrorException check_parameters(testeq, NamedTuple())
@test_throws ErrorException check_parameters(testeq, 1)

@test_throws ErrorException datatype(testeq, NamedTuple())
@test_throws ErrorException arrtype(testeq, NamedTuple())

@test hassolution(testeq) == false
@test hasvectorfield(testeq) == false
@test hasprimary(testeq) == false
@test hassecondary(testeq) == false

@test hasinvariants(testeq) == false
@test hasparameters(testeq) == false
@test hasperiodicity(testeq) == false

@test hashamiltonian(testeq) == false
@test haslagrangian(testeq) == false
