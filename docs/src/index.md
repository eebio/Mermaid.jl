# Mermaid.jl

```@meta
CurrentModule = Mermaid
```

## Summary / README

```@docs
Mermaid
```

## Installation (Coming soon)

Mermaid can be installed from Julia with:

```julia
using Pkg; Pkg.add("Mermaid")
```

## FAQ

**Q: How can I use Mermaid to run a hybrid model simulation?**

**A:** Check out the [Tutorial](@ref), it covers how to set up and run a hybrid simualtion between an Agent-based model and an ODE system.

**Q: What if Mermaid doesn't have a component for the type of model I want to use?**

**A:** Check out the [Mermaid Interface](@ref), you can easily define new components for many different simulation tools.

**Q: Where can I view more examples?**

**A:** In the Examples section of course. We have examples showing all of the more advanced features of Mermaid including:

* duplicated components,
* out of sync computation,
* external (non-Julia) components,
* surrogate model approximations,
* mapping between spatial components on different resolutions,
* [integration with ModelingToolkit.jl](@ref "ModelingToolkit").

**Q: Is Mermaid the right tool for me?**
**A:** That depends on what type of hybrid simulation you are looking to run. Luckily, you can check out [Is Mermaid right for me?](@ref) and find out.
