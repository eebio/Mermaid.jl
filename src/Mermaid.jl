module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

# Imports, Usings and Reexports
import CommonSolve: solve!, solve, init, step!
using OrderedCollections: OrderedDict
export OrderedDict

# Exports
export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export AbstractComponentIntegrator
export DEComponent, DuplicatedComponent, MOLComponent, AgentsComponent, SurrogateComponent
export Connector, ConnectedVariable, MermaidProblem, MermaidIntegrator
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve, init, step!
export getstate, setstate!, gettime, settime!
export name, time_step, variables
export fullname, runconnection, runconnection!

# Include src files
include("abstracts.jl")
include("connections.jl")
include("solutions.jl")
include("mermaidProblems.jl")
include("solvers.jl")
include("Duplicated.jl")
include("extensions.jl")

# Documentation
"""
    setstate!(comp::AbstractComponentIntegrator, state)
    setstate!(comp::AbstractComponentIntegrator, key, value)

Set the state of a component.

# Arguments
- `comp::AbstractComponentIntegrator`: The component whose state is to be set.
- `state`: The new state to set for the entire component.
- `key`: The key specifying which part of the component's state to set.
- `value`: The value to set for the specified part of the component's state.
"""
function setstate! end

"""
    getstate(comp::AbstractComponentIntegrator)
    getstate(comp::AbstractComponentIntegrator, key)

Retrieve the state of a component.

# Arguments
- `comp::AbstractComponentIntegrator`: The component whose state is to be retrieved.
- `key`: The key specifying which part of the component's state to retrieve.

"""
function getstate end

"""
    variables(comp::AbstractComponentIntegrator)
    variables(comp::AbstractComponent)

Retrieve the variable names of a component.

# Arguments
- `comp::Union{AbstractComponent, AbstractComponentIntegrator}`: The component (or component
    integrator) whose variable names are to be retrieved.

# Returns
- A collection of variable names (as strings) associated with the component. This includes
    all special variables such as `#time` and `#model` if applicable.
"""
function variables end

"""
    step!(int::AbstractMermaidIntegrator)
    step!(int::AbstractComponentIntegrator)

Advance the state of the integrator `int` by one time step.

# Arguments
- `int::Union{AbstractMermaidIntegrator, AbstractComponentIntegrator}`: The integrator to
    advance.
"""
function step! end

"""
    init(prob::MermaidProblem, alg::MermaidSolver; save_vars=[])
    init(comp::AbstractComponent)

Initialises an integrator ([MermaidIntegrator](@ref) or [AbstractComponentIntegrator](@ref))
    for the given [MermaidProblem](@ref)/[AbstractComponent](@ref).

# Arguments
- `prob::MermaidProblem`: The hybrid problem to be solved.
- `alg::AbstractMermaidSolver`: The [AbstractMermaidSolver](@ref) algorithm to be used for
    solving the problem.
- `comp::AbstractComponent`: The component to be initialised.

# Keyword Arguments
- `save_vars=[]`: Which variables to save in the solution. Defaults to all non-special
    variables.

# Returns
- `MermaidIntegrator`: The initialized integrator for the problem.
"""
function init end

end
