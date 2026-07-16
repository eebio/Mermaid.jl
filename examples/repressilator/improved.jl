module Improved
    using Catalyst
    using StochasticDiffEq
    import Distributions: Geometric
    using Polynomials

    # Constants
    λ = 60
    K = [5, 10, 10]
    n = 2
    b̄ = 10
    N̄ₒ = 10
    N̄ₜ = 40

    @register_symbolic Geometric(b)
    @parameters b
    m = rand(Geometric(1/b)) + 1

    function root(Pt, Npro, K)
        p = Polynomial([-K^1.5*Pt, 0, K^1.5, 2Npro-Pt, 0, 1])
        f = x -> isreal(x) && 0<=real(x)<=Pt
        r = roots(p)
        r = filter(f, r)
        if isempty(r)
            if Pt < 1e-3
                return 0.0
            else
                return Pt
            end
        end
        return real(first(r))^2
    end

    @register_symbolic root(Pt, Npro, K)

    # Nₒ and Nₜ are the numbers of repressilator and titration sponge plasmids (means of N̄ₒ and N̄ₜ)
    # P₁, P₂, P₃ are the three repressors (cI, lacI, tetR)
    # gfp is the green fluorescent protein reporter
    function sde_improved_repressilator()
        r = @reaction_network begin
            @parameters n K[1:3]
            @species gr(t) V(t)
            a*gr * N̄ₒ*V, ∅ --> Nₒ
            a*gr * N̄ₜ*V, ∅ --> Nₜ
            a*gr * (hillr(root(P₃/V, (Nₒ + Nₜ)/V, K[3]), λ * (Nₒ + Nₜ)/V, K[3], n))*V, ∅ --> b̄ * P₁
            a*gr * (hillr(root(P₁/V, Nₒ/V, K[1]), λ * Nₒ/V, K[1], n))*V, ∅ --> b̄ * P₂
            a*gr * (hillr(root(P₂/V, Nₒ/V, K[2]), λ * Nₒ/V, K[2], n))*V, ∅ --> b̄ * P₃
            a*gr * (hillr(root(P₃/V, (Nₒ + Nₜ)/V, K[3]), λ * (Nₒ + Nₜ)/V, K[3], n))*V, ∅ --> b̄ * gfp
            b*gr, Nₒ --> ∅ # Degradation is only mimicing dilution
            b*gr, Nₜ --> ∅
            b*gr, P₁ --> ∅
            b*gr, P₂ --> ∅
            b*gr, P₃ --> ∅
            b*gr, gfp --> ∅
        end

        u0 = [:P₁ => 0.0, :P₂ => 0.0, :Nₒ => 10.0, :Nₜ => 40.0, :P₃ => 1000.0, :gfp => 0.0, :gr => 1.0, :V => 1.0]
        tspan = (0., 30.)
        ps = [:K => K, :N̄ₒ => N̄ₒ, :N̄ₜ => N̄ₜ, :λ => λ, :n => n, :b̄ => b̄, :a => 15.0, :b => 4.0]
        sde = SDEProblem(r, u0, tspan, ps; structural_simplify = true)
        return sde
    end
end
