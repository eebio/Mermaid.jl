module Mermaid
using CommonSolve
using DifferentialEquations
using OrdinaryDiffEqCore

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent, PDEComponent, MermaidProblem
export Connector
export AbstractMermaidSolver, MinimumTimeStepper
export solve!, solve
include("types.jl")
include("functions.jl")
include("DiffEq.jl")
include("PDE.jl")
include("solvers.jl")

end
