using CommonSolve

"""
    MermaidProblem <: AbstractMermaidProblem
    MermaidProblem(;
        components::Vector{AbstractComponent},
        connectors::Vector{AbstractConnector},
        max_t::Float64=1.0,
        timescales::Vector{Float64}=ones(length(components)))

Defines a Mermaid hybrid simulation problem.

# Keyword Arguments
- `components::Vector{<:AbstractComponent}`: Vector of Components.
- `connectors::Vector{<:AbstractConnector}`: Vector of [Connector](@ref).
- `max_t::Float64`=1.0: Maximum simulation time. Defaults to 1.0.
- `timescales::Vector{Float64}`=ones(length(components)): Timescales for each component.
    Each component's timescale will be multiplied by the component's time to convert it to
    the universal simulation time.
"""
@kwdef struct MermaidProblem <: AbstractMermaidProblem
    components::Vector{AbstractComponent}
    connectors::Vector{AbstractConnector}
    max_t::Float64 = 1.0
    timescales::Vector{Float64} = ones(length(components))
end

"""
    MermaidIntegrator <: AbstractMermaidIntegrator
    MermaidIntegrator(
        integrators::Vector{AbstractComponentIntegrator},
        connectors::Vector{AbstractConnector},
        maxt::Float64,
        currtime::Float64,
        alg::AbstractMermaidSolver,
        save_vars::Vector{String})

Defines the integrator for a Mermaid hybrid simulation.

# Fields
- `integrators::Vector{<:AbstractComponentIntegrator}`: Vector of component integrators.
- `connectors::Vector{<:AbstractConnector}`: Vector of [Connector](@ref) between components.
- `maxt::Float64`: Maximum simulation time.
- `currtime::Float64`: Current simulation time.
- `alg::AbstractMermaidSolver`: The Mermaid solver algorithm to be used.
- `save_vars::Vector{String}`: Variables to be saved during the simulation.
"""
mutable struct MermaidIntegrator{X <: AbstractMermaidSolver} <: AbstractMermaidIntegrator
    integrators::Vector{<:AbstractComponentIntegrator}
    connectors::Vector{<:AbstractConnector}
    maxt::Float64
    currtime::Float64
    alg::X
    save_vars::Vector{<:AbstractString}
    timescales::Vector{Float64}
end

function CommonSolve.init(
        prob::AbstractMermaidProblem, alg::AbstractMermaidSolver; save_vars = String[])
    # Initialize the solver
    integrators = [init(c) for c in prob.components]
    return MermaidIntegrator(
        integrators, prob.connectors, prob.max_t, 0.0, alg, save_vars, prob.timescales)
end

function CommonSolve.step!(merInt::AbstractMermaidIntegrator)
    step!(merInt, merInt.alg)
end

"""
    solve!(merInt::AbstractMermaidIntegrator)

Solves the problem using the [MermaidIntegrator](@ref). This handles all the message passing
    and calls step! on the [MermaidIntegrator](@ref).

# Arguments
- `merInt::AbstractMermaidIntegrator`: The integrator to be solved.

# Returns
- `MermaidSolution`: The [solution](@ref MermaidSolution) of the problem.
"""
function CommonSolve.solve!(merInt::AbstractMermaidIntegrator)
    sol = MermaidSolution(merInt)
    update_solution!(sol, merInt)
    while merInt.currtime < merInt.maxt
        step!(merInt)
        update_solution!(sol, merInt)
    end
    return sol
end

function getstate(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable)
    # Get the state of the component based on the key
    for integrator in merInt.integrators
        if name(integrator) == key.component
            return getstate(integrator, key)
        end
    end
end

function gettime(merInt::AbstractMermaidIntegrator)
    # Get the current time of the integrator
    return merInt.currtime
end

function setstate!(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable, value)
    # Set the state of the component based on the key
    for integrator in merInt.integrators
        if name(integrator) == key.component
            setstate!(integrator, key, value)
            return nothing
        end
    end
end

"""
    gettime(merInt::AbstractComponentIntegrator)

Get the current time of the integrator.

# Arguments
- `int::AbstractComponentIntegrator`: The integrator whose time is to be retrieved.

# Returns
- The current time of the integrator.
"""
function gettime(int::AbstractComponentIntegrator)
    getstate(
        int, ConnectedVariable(name(int), "#time", nothing, nothing))
end

"""
    settime!(merInt::AbstractComponentIntegrator, t)

Set the current time of the integrator.

# Arguments
- `int::AbstractComponentIntegrator`: The integrator whose time is to be set.
- `t`: The time to set.
"""
function settime!(int::AbstractComponentIntegrator, t)
    setstate!(int,
        ConnectedVariable(name(int), "#time", nothing, nothing), t)
end

"""
    variables(integrator::AbstractComponent)
    variables(integrator::AbstractComponentIntegrator)

Get the variables names of the component or integrator that can be accessed through getstate
    and setstate!. This includes special variables like `#time`.
"""
variables(integrator::AbstractComponentIntegrator) = variables(integrator.component)

function getstate(int::AbstractComponentIntegrator; copy = false)
    if copy
        return deepcopy(getstate(int))
    else
        return getstate(int)
    end
end

"""
    time_step(int::AbstractComponent)
    time_step(comp::AbstractComponentIntegrator)

Get the proposed time step of the integrator or component. It can depend on the current
    state.
"""
time_step(comp::AbstractComponent) = comp.time_step
time_step(int::AbstractComponentIntegrator) = time_step(int.component)

"""
    name(int::AbstractComponentIntegrator)
    name(comp::AbstractComponent)

Get the name of the integrator or component.
"""
name(comp::AbstractComponent) = comp.name
name(int::AbstractComponentIntegrator) = name(int.component)
