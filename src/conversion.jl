
macro define(name, definition)
    quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


@define _create_pode_argument_views begin
    n = div(length(eachindex(x)), 2)
    q = @view x[eachindex(x)[  1:n ]]
    p = @view x[eachindex(x)[n+1:2n]]
    q̇ = @view ẋ[eachindex(ẋ)[  1:n ]]
    ṗ = @view ẋ[eachindex(ẋ)[n+1:2n]]
end


extend_periodicity(equ::AbstractEquationPODE) = periodicity(equ) == NullPeriodicity() ? periodicity(equ) : vcat(periodicity(equ), zero(periodicity(equ)))

convert_periodicity(::Union{Type{ODE}, Type{SODE}}, equ::Union{PODE, HODE}) = extend_periodicity(equ)
convert_periodicity(::Union{Type{ODE}, Type{SODE}}, prob::Union{PODEProblem, HODEProblem}) = extend_periodicity(equation(prob))


function Base.convert(::Type{ODEProblem}, prob::Union{PODEProblem{DT,TT,AT}, HODEProblem{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = vcat(prob.ics.q, prob.ics.p)

    # extend periodicity
    ode_periodicity = convert_periodicity(ODE, prob)

    v = (ẋ, t, x, params) -> begin
        @_create_pode_argument_views
        equation(prob).v(q̇, t, q, p, params)
        equation(prob).f(ṗ, t, q, p, params)
    end

    ODEProblem(v, prob.tspan, prob.tstep, x₀; parameters=prob.parameters, periodicity=ode_periodicity)
    # TODO: Convert invariants and pass to ODE
    # TODO: For HODE append (h=equ.h,) to invariants
end

function Base.convert(::Type{SODEProblem}, prob::Union{PODEProblem{DT,TT,AT}, HODEProblem{DT,TT,AT}}) where {DT, TT, AT <: AbstractVector}
    # concatenate initial conditions
    x₀ = vcat(prob.ics.q, prob.ics.p)

    # extend periodicity
    ode_periodicity = convert_periodicity(SODE, prob)

    v₁ = (ẋ, t, x, params) -> begin
        @_create_pode_argument_views
        equation(prob).v(q̇, t, q, p, params)
    end
    v₂ = (ẋ, t, x, params) -> begin
        @_create_pode_argument_views
        equation(prob).f(ṗ, t, q, p, params)
    end

    SODEProblem((v₁, v₂), prob.tspan, prob.tstep, x₀; parameters=parameters(prob), periodicity=ode_periodicity)
    # TODO: Convert invariants and pass to SODE
end

function Base.convert(::Type{PODEProblem}, prob::HODEProblem)
    PODEProblem(equation(prob).v, equation(prob).f, prob.tspan, prob.tstep, prob.ics.q, prob.ics.p;
                invariants=invariants(equation(prob)), parameters=parameters(prob), periodicity=periodicity(equation(prob)))
    # TODO: Append (h=equ.h,) to invariants
end

function Base.convert(::Type{IODEProblem}, prob::LODEProblem)
    IODEProblem(equation(prob).ϑ, equation(prob).f, equation(prob).g, prob.tspan, prob.tspep, prob.ics.q, prob.ics.p, prob.ics.λ;
                v̄=equation(prob).v̄, f̄=equation(prob).f̄, invariants=invariants(equation(prob)), parameters=parameters(prob), periodicity=periodicity(equation(prob)))
end


function get_invariants(equ::Union{ODE,SODE,DAE})
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t,q) -> inv(t, q, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end

function get_invariants(equ::Union{IODE,LODE,IDAE,LDAE})
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t, q, v) -> inv(t, q, v, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end

function get_invariants(equ::Union{PODE,HODE,PDAE,PDAE})#,SPDAE
    if hasinvariants(equ)
        keys = ()
        invs = ()
        for (key, inv) in pairs(equ.invariants)
            keys = (keys..., key)
            invs = (invs..., hasparameters(equ) ? (t, q, p) -> inv(t, q, p, equ.parameters) : inv)
        end
        return NamedTuple{keys}(invs)
    else
        return NamedTuple()
    end
end
