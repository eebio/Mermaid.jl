module Mermaid
using CommonSolve
using DifferentialEquations
using OrdinaryDiffEqCore

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent, MermaidProblem
export Connector
export AbstractMermaidSolver, MinimumTimeStepper
export solve!, solve
include("types.jl")
include("functions.jl")
include("DiffEq.jl")
include("solvers.jl")

end
