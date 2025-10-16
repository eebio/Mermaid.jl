module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

# Imports, Usings and Reexports
using CommonSolve
using OrderedCollections: OrderedDict
export OrderedDict

# Exports
export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ComponentIntegrator
export DEComponent, DuplicatedComponent, PDEComponent, AgentsComponent, SurrogateComponent
export Connector, ConnectedVariable, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve, init
export getstate, setstate!, gettime, settime!, step!, variables
export inputsandoutputs, parsevariable, fullname

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

end
