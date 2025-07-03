# ModelingToolkit Integration

Similarly to Mermaid, [ModelingToolkit](@extref ModelingToolkit index) also allows specifying models as components and then connecting them together.
However, the methods used are different.
For an overview of these differences, you can see [Is Mermaid right for me?](@ref).
In summary, if it is possible to connect your components through ModelingToolkit, it is likely better to do that rather than through Mermaid.
In order to simplify connected ModelingToolkit systems in Mermaid, you can create `ModelingToolkitComponents` for all of your ModelingToolkit models, and specify the connections through Mermaid.
Mermaid will then figure out which `ModelingToolkitComponents` can be connected together in the [MermaidProblem](@ref), perform the connections, generate the model and solve it.

## ModelingToolkitComponents
TODO
