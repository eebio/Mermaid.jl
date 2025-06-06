module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

using CommonSolve
using DifferentialEquations
using OrdinaryDiffEqCore

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent, PDEComponent, AgentsComponent
export DuplicatedComponent
export Connector, ConnectedVariable, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve
include("types.jl")
include("functions.jl")
include("DiffEq.jl")
include("PDE.jl")
include("Duplicated.jl")
include("Agents.jl")
include("solvers.jl")

end
