module Mermaid
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) Mermaid

using CommonSolve
using OrderedCollections: OrderedDict
export OrderedDict

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ComponentIntegrator
export DEComponent, DuplicatedComponent, PDEComponent, AgentsComponent, SurrogateComponent
export Connector, ConnectedVariable, MermaidProblem
export AbstractMermaidSolver, MinimumTimeStepper
export MermaidSolution
export solve!, solve, init
export getstate, setstate!, gettime, settime!, step!, variables
export inputsandoutputs, parsevariable, fullname

include("abstracts.jl")
include("connections.jl")
include("solutions.jl")
include("mermaidProblems.jl")
include("solvers.jl")
include("Duplicated.jl")
include("extensions.jl")

end
