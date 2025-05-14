module Mermaid
using CommonSolve
using DifferentialEquations
using OrdinaryDiffEqCore

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent, PDEComponent, AgentsComponent
export Connector, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export solve!, solve
include("types.jl")
include("functions.jl")
include("DiffEq.jl")
include("PDE.jl")
include("Agents.jl")
include("solvers.jl")

end
