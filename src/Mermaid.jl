module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

using CommonSolve
using OrderedCollections

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export DEComponent, PDEComponent, AgentsComponent, SurrogateComponent
export DuplicatedComponent
export Connector, ConnectedVariable, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve
export inputsandoutputs
include("types.jl")
include("functions.jl")
include("DiffEq.jl")
include("PDE.jl")
include("Duplicated.jl")
include("Agents.jl")
include("Surrogates.jl")
include("solvers.jl")

end
