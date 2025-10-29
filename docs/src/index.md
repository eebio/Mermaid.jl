# Mermaid.jl

```@meta
CurrentModule = Mermaid
```

## Summary / README

```@docs
Mermaid
```

## Installation

Mermaid can be installed from Julia with:

```julia
using Pkg; Pkg.add("Mermaid")
```

## FAQ

**Q: How can I use Mermaid to run a hybrid model simulation?**

**A:** Check out the [Tutorial](@ref), it covers how to set up and run a hybrid simualtion between an Agent-based model and an ODE system.

**Q: What if Mermaid doesn't have a component for the type of model I want to use?**

**A:** Check out the [Mermaid Interface](@ref), you can easily define new components for many different simulation tools.

**Q: The Tutorial didn't cover what I wanted. Where can I view more examples?**

**A:** In the Examples section of course. We have examples showing all of the more advanced features of Mermaid including:

* [duplicated components](@ref "Advanced Duplicated Components"),
* [out of sync computation](@ref "Out of sync computation"),
* [external (non-Julia) components](@ref "External Components"),
* [surrogate model approximations](@ref "Surrogates"),
* [mapping between spatial components on different resolutions](@ref "Spatial maps"),
* [integration with ModelingToolkit.jl](@ref "ModelingToolkit Integration").

**Q: Is Mermaid the right tool for me?**

**A:** That depends on what type of hybrid simulation you are looking to run. Luckily, you can check out [Is Mermaid right for me?](@ref) and find out.
