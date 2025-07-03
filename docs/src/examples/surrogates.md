# Surrogates

Sometimes, a model can be very expensive to solve.
In other instances, we may be required to solve the same model many times, such as with [Duplicated Components](@ref "Advanced Duplicated Components").
In these instances, it can be useful to generate a surrogate machine learning model which approximates the output, and can do it faster than model solving.

Surrogates of components can be generated in Mermaid.
The component will train a surrogate to learn the `step!` function.
