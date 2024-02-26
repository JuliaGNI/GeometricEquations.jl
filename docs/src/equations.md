```@meta
CurrentModule = GeometricEquations
```

# Equations

Equations hold a number of functions determining vector fields, constraints, invariants, the Hamiltonian or Lagrangian and the symplectic structure.

There exists a hierarchy of abstract data types. At the top stands the `GeometricEquation` type. From that several types derive for ODEs, DAEs, SDEs and their partitioned counterparts:

- `AbstractEquationODE`
- `AbstractEquationDAE`
- `AbstractEquationSDE`
- `AbstractEquationPODE`
- `AbstractEquationPDAE`
- `AbstractEquationPSDE`

Concrete implementations of equation types should be subtypes of one of these types or abstract subtypes thereof.


Equations have several traits:

| Function                             | Description   |
|:------------------------------------ |:------------- |
| `hassolution(::GeometricEquation)`   |               | 
| `hasvectorfield(::GeometricEquation)`|               |   
| `hasprimary(::GeometricEquation)`    |               | 
| `hassecondary(::GeometricEquation)`  |               | 
| `hasinvariants(::GeometricEquation)` |               |  
| `hasparameters(::GeometricEquation)` |               |  
| `hasperiodicity(::GeometricEquation)`|               |   
| `hashamiltonian(::GeometricEquation)`|               |   
| `haslagrangian(::GeometricEquation)` |               |  

Equations should implement some convenience functions for accessing their data:

| Function                             | Description   |
|:------------------------------------ |:------------- |
| `functions(::GeometricEquation)`     |               |
| `solutions(::GeometricEquation)`     |               |
| `invariants(::GeometricEquation)`    |               |
| `parameters(::GeometricEquation)`    |               |
| `periodicity(::GeometricEquation)`   |               |

For each equation type, there should be some methods for functions implemented, that are used to check the validity of initial conditions or the functions provided for the vector fields, etc.:

`check_initial_conditions(equ::GeometricEquation, ics::NamedTuple)`
`check_parameters(equ::GeometricEquation, params::NamedTuple)`



## Ordinary Differential Equations

```@docs
ODE
PODE
HODE
IODE
LODE
SODE
```


## Differential Algebraic Equations

```@docs
DAE
PDAE
HDAE
IDAE
LDAE
```


## Stochastic Differential Equations

```@docs
SDE
PSDE
SPSDE
```
