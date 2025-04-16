module Mermaid
using CommonSolve
using DifferentialEquations
using OrdinaryDiffEqCore

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent, MermaidProblem
export MermaidSolver
export solve!, solve
include("types.jl")
include("functions.jl")

end
