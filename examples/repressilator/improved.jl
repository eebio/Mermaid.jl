using Catalyst
using DifferentialEquations, StochasticDiffEq
using Plots
using JumpProcesses
import Distributions: Geometric
using Polynomials
using Random

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

function jump_improved_repressilator()
    r = @reaction_network begin
        @species gr(t)
        @parameters K[1:3] P₁_free P₂_free P₃_free
        gr * N̄ₒ, ∅ --> Nₒ
        gr * N̄ₜ, ∅ --> Nₜ
        gr * (hillr(root(P₃, Nₒ+Nₜ, K[3]), λ * (Nₒ+Nₜ), K[3], n)), ∅ --> $m * P₁
        gr * (hillr(root(P₁, Nₒ, K[1]), λ * Nₒ, K[1], n)), ∅ --> $m * P₂
        gr * (hillr(root(P₂, Nₒ, K[2]), λ * Nₒ, K[2], n)), ∅ --> $m * P₃
        gr * (hillr(root(P₃, Nₒ+Nₜ, K[3]), λ * (Nₒ * Nₜ), K[3], n)), ∅ --> $m * gfp
        gr, Nₒ --> ∅ # Degradation is only mimicing dilution
        gr, Nₜ --> ∅
        gr, P₁ --> ∅
        gr, P₂ --> ∅
        gr, P₃ --> ∅
        gr, gfp --> ∅
    end

    u0 = [:P₁ => 160.0, :P₂ => 55.0, :Nₒ => 10.0, :Nₜ => 40.0, :P₃ => 2400.0, :gfp => 0.0, :gr => 1.0]
    tspan = (0., 60.)
    ps = [:K => K, :N̄ₒ => N̄ₒ, :N̄ₜ => N̄ₜ, :λ => λ, :n => n, :b => b̄, :P₁_free => 50.0, :P₂_free => 40.0, :P₃_free => 30.0]
    jsys = convert(JumpSystem, r; combinatoric_ratelaws=false)
    jsys = complete(jsys)
    dprob = DiscreteProblem(jsys, u0, tspan, ps)
    jprob = JumpProblem(jsys, dprob, Direct())
    return jprob, r
end

function plot_jump()
    jprob, r = jump_improved_repressilator()
    sol = solve(jprob, SSAStepper())
    display(Plots.plot(sol, vars=[:P₁, :P₂, :P₃], xlabel="Time", ylabel="Concentration", title="Repressilator Protein Dynamics (Gillespie)", linewidth=2))
    display(Plots.plot(sol, vars=[:Nₒ, :Nₜ], xlabel="Time", ylabel="Concentration", title="Plasmid Concentrations (Gillespie)", linewidth=2))
    Plots.plot(sol, vars=[:gfp], xlabel="Time", ylabel="Concentration", title="Repressilator Protein Dynamics (Gillespie)", linewidth=2)

    sol = solve(jprob, SSAStepper())
    display(Plots.plot!(sol, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))
end

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

function plot_sde()
    sde = sde_improved_repressilator()
    method = EM()
    sol_sde = solve(sde, method, dt=0.01)

    display(Plots.plot(sol_sde, vars=[:P₁, :P₂, :P₃], xlabel="Time", ylabel="Concentration", title="Repressilator Protein Dynamics (SDE)", linewidth=2))
    display(Plots.plot(sol_sde, vars=[:Nₒ, :Nₜ], xlabel="Time", ylabel="Concentration", title="Plasmid Concentrations (SDE)", linewidth=2))
    Plots.plot(sol_sde, vars=[:gfp], xlabel="Time", ylabel="Concentration", title="Repressilator Protein Dynamics (SDE)", linewidth=2)

    sol_sde = solve(sde, method, dt=0.01)
    display(Plots.plot!(sol_sde, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))
end

function final_plot()
    Random.seed!(0) # The seeding here is crazy, different between Shift+Enter and running the lines individually, not running the same way twice, etc., there is a seed option in StochasticDiffEq but doesn't seem to work
    prob = sde_improved_repressilator()
    Random.seed!(0)
    sol_sde = solve(prob, EM(), dt=0.1, save_idxs=[:gfp, :P₁, :P₂, :P₃], saveat=0.1)
    display(sol_sde[end])
    average_gfp = sol_sde[:gfp]
    for i in 1:30
        Random.seed!(i)
        sol2 = solve(prob, EM(), dt=0.1, save_idxs=[:gfp], saveat=0.1)
        average_gfp .+= sol2[:gfp]
    end
    average_gfp ./= 31
    p1 = Plots.plot(sol_sde.t, sol_sde[:gfp], label="GFP", xlabel="Time", ylabel="Concentration", title="Improved Repressilator GFP", linewidth=2, ylims=(0, 1.3 * maximum(sol_sde[:gfp])))
    Plots.plot!(p1, sol_sde.t, average_gfp, linewidth=2, linestyle=:dash, label="Average GFP")

    p2 = Plots.plot(sol_sde.t, sol_sde[:P₁], label="tetR", xlabel="Time", ylabel="Concentration", title="Proteins", linewidth=2)
    Plots.plot!(p2, sol_sde.t, sol_sde[:P₂], label="cI", linewidth=2)
    Plots.plot!(p2, sol_sde.t, sol_sde[:P₃], label="lacI", linewidth=2)

    # Plots.plot!(p1, xticks=0:250:750)
    # Plots.plot!(p2, xticks=0:250:750)

    l = @layout [a b]
    display(Plots.plot(p1, p2, layout=l, size=(750, 350), margin=2*Plots.mm))
    savefig("repressilator_plots_improved.png")
end

#plot_jump()
#plot_sde()
#final_plot()
