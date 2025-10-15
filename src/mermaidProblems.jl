using CommonSolve

"""
    MermaidProblem

Struct for defining a Mermaid problem.
This struct contains the components, connectors and other properties of the hybrid simulation.

# Fields
- `components::Vector`: Vector of Components.
- `connectors::Vector{Connector}`: Vector of [Connector](@ref).
- `max_t::Float64`=1.0: Maximum simulation time.
"""
@kwdef struct MermaidProblem <: AbstractMermaidProblem
    components::Vector{AbstractComponent}
    connectors::Vector{AbstractConnector}
    max_t::Float64 = 1.0
end

"""
    MermaidIntegrator

Mutable struct for integrating a hybrid [MermaidProblem](@ref).
This struct holds all the [ComponentIntegrators](@ref ComponentIntegrator) and [Connectors](@ref Connector) to store the current state of the hybrid simulation.

# Fields
- `integrators::Vector`: Vector of [ComponentIntegrator](@ref).
- `connectors::Vector{Connector}`: Vector of [Connector](@ref).
- `maxt::Float64`: Maximum simulation time.
- `currtime::Float64`: Current simulation time.
- `alg::AbstractMermaidSolver`: Algorithm used for integration.
- `save_vars::Vector{String}`: Variables to save during integration.
"""
mutable struct MermaidIntegrator <: AbstractMermaidIntegrator
    integrators::Vector{AbstractComponentIntegrator}
    connectors::Vector{AbstractConnector}
    maxt::Float64
    currtime::Float64
    alg::AbstractMermaidSolver
    save_vars::Vector{String}
end

"""
    init(prob::MermaidProblem, alg::MermaidSolver; kwargs...)

Creates a [MermaidIntegrator](@ref) from a [MermaidProblem](@ref).

# Arguments
- `prob::MermaidProblem`: The hybrid problem to be solved.
- `alg::AbstractMermaidSolver`: The [AbstractMermaidSolver](@ref) algorithm to be used for solving the problem.
- `kwargs...`: Additional keyword arguments for the solver.

# Returns
- `MermaidIntegrator`: The initialized integrator for the problem.
"""
function CommonSolve.init(prob::AbstractMermaidProblem, alg::AbstractMermaidSolver; save_vars=[])
    # Sort connectors so that the ones with # in outputs are first
    connectors = sort(prob.connectors; by=x -> any([contains(conn.variable, "#") for conn in x.outputs]), rev=true)
    # Initialize the solver
    integrators = Vector{Any}()
    for c in prob.components
        integrator = CommonSolve.init(c, connectors)
        push!(integrators, integrator)
    end
    return MermaidIntegrator(integrators, connectors, prob.max_t, 0.0, alg, save_vars)
end

"""
    step!(int::MermaidIntegrator, dt::Float64)

Steps the integrator for the given time step. The method is defined by the `alg` field of the [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to be stepped.
- `dt::Float64`: The time step for the integrator.
"""
function CommonSolve.step!(merInt::AbstractMermaidIntegrator, dt)
    merInt.alg(merInt, dt)
end

"""
    solve!(int::MermaidIntegrator)

Solves the problem using the [MermaidIntegrator](@ref).
This handles all the message passing and calls step! on the [MermaidIntegrator](@ref).

# Arguments
- `int::MermaidIntegrator`: The integrator to be solved.

# Returns
- `MermaidSolution`: The [solution](@ref MermaidSolution) of the problem.
"""
function CommonSolve.solve!(int::AbstractMermaidIntegrator)
    t = [0.0]
    dt = minimum([i.component.time_step for i in int.integrators]) # Minimum isnt sufficient to guarantee we don't jump over anything
    sol = MermaidSolution(int)
    update_solution!(sol, int)
    while int.currtime < int.maxt
        CommonSolve.step!(int, dt)
        update_solution!(sol, int)
    end
    return sol
end

"""
    getstate(merInt::MermaidIntegrator, key::ConnectedVariable)

Retrieve the state of a component within a [MermaidIntegrator](@ref) based on the provided [ConnectedVariable](@ref).

# Arguments
- `merInt::MermaidIntegrator`: The integrator containing multiple component integrators.
- `key::ConnectedVariable`: The key specifying which component's state to retrieve.

# Returns
- The state associated with the specified key, or `nothing` if not found.
"""
function getstate(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable)
    # Get the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
            return getstate(integrator, key)
        end
    end
end

"""
    setstate!(merInt::MermaidIntegrator, key::ConnectedVariable, value)

Set the state of a component integrator inside a [MermaidIntegrator](@ref) based on the provided key and value.

# Arguments
- `merInt::MermaidIntegrator`: The integrator containing multiple component integrators.
- `key::ConnectedVariable`: The key specifying which component's state to set.
- `value`: The value to assign to the specified state.
"""
function setstate!(merInt::AbstractMermaidIntegrator, key::AbstractConnectedVariable, value)
    # Set the state of the component based on the key
    for integrator in merInt.integrators
        if integrator.component.name == key.component
            setstate!(integrator, key, value)
            return nothing
        end
    end
end

gettime(integrator::AbstractComponentIntegrator) = getstate(integrator, ConnectedVariable(integrator.component.name, "#time", nothing, nothing))
settime!(integrator::AbstractComponentIntegrator, t) = setstate!(integrator, ConnectedVariable(integrator.component.name, "#time", nothing, nothing), t)

variables(integrator::AbstractComponentIntegrator) = variables(integrator.component)
