module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

using CommonSolve
using OrderedCollections

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ComponentIntegrator
export DEComponent, DuplicatedComponent, PDEComponent, AgentsComponent, SurrogateComponent
export Connector, ConnectedVariable, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve, init
export getstate, setstate!, gettime, settime!, step!, variables
export inputsandoutputs, parsevariable, fullname
include("types.jl")
include("functions.jl")
include("Duplicated.jl")
include("solvers.jl")
include("extensions.jl")

end
