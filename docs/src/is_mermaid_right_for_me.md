# Is Mermaid right for me?

This page will walk through the benefits and drawbacks of using Mermaid, what it can do, what it does well, what it can't do and what it does poorly.

## Comparison against ModelingToolkit.jl

ModelingToolkit.jl is an acausal modeling language, where you specify symbolic forms of your models, which can then be connected together, and compiled into a DifferentialEquations.jl problem.

In short, if your whole system can be expressed in ModelingToolkit and solved (*mostly* equivalently, the combined system can be solved in DifferentialEquations or JumpProcesses), then it is probably better solved through ModelingToolkit rather than Mermaid. The acausal connections lead to better accuracy and stability of the final solve.

The most notable exception to this is multirate solvers. Mermaid naturally has components solved through different solvers, on different (possibly adaptive) timescales. Multirate and multisolver functionality is limited in DifferentialEquations.

!!! warning "Multi-rate solvers"
    Using multirate solver for differential equations typically comes with a lower order accuracy on the connected variables. If you are unsure if a multirate solver is suitable, it is likely better to favour the solver algorithms in DifferentialEquations and use ModelingToolkit to perform the connections.

If only part of the whole system can be expressed in ModelingToolkit, you may still be able to benefit from ModelingToolkit features. Combining the relevant components together through ModelingToolkit and then converting that `System` to a component to connect through Mermaid allows you to benefit from all ModelingToolkit features possible, while still having the flexibility to connect outside the ModelingToolkit ecosystem.
