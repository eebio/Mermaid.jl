module Mermaid
using DifferentialEquations

export AbstractComponent, AbstractTimeDependentComponent, AbstractTimeIndependentComponent
export ODEComponent
export solve!
include("types.jl")
include("functions.jl")

end
