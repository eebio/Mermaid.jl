using Catalyst
using DifferentialEquations, StochasticDiffEq
using Plots

# Constants
s = 10.0^4
γ_max = 1260.0
νₜ_max = 726.0
νₘ_max = 5800.0
wᵣ_max = 930.0
wₕ_max = 949.0
θₜ = 4.38
θₘ = 4.38
θₕ = 4.38
θᵣ = 427.0
k = 0.0095
Mₛ = 10.0^8
Kᵧ = 7.0
Kₜ = 1000.0
Kₘ = 1000.0
Kₕ = 152219.0
wₘ_max = 4.14
wₜ_max = 4.14
dₘ = 0.1
hₕ = 4.0
nᵣ = 7459
nₕ = 300
nₜ = 300
nₘ = 300
kᵤ = 1.0
nₛ = 1

function get_growth_model()

    growth = @reaction_network begin
        @species s
        @parameters s γ_max νₜ_max νₘ_max wᵣ_max wₕ_max θₜ θₘ θₕ θᵣ k Mₛ Kᵧ Kₜ Kₘ Kₕ wₘ_max wₜ_max dₘ hₕ nᵣ nₕ nₜ nₘ kᵤ nₛ
        # Mass action reactions
        (wᵣ, λ + dₘ), ∅ <--> mᵣ
        (wₜ, λ + dₘ), ∅ <--> mₜ
        (wₘ, λ + dₘ), ∅ <--> mₘ
        (wₕ, λ + dₘ), ∅ <--> mₕ
        (k, kᵤ), pᵣ + mᵣ <--> cᵣ
        (k, kᵤ), pᵣ + mₜ <--> cₜ
        (k, kᵤ), pᵣ + mₘ <--> cₘ
        (k, kᵤ), pᵣ + mₕ <--> cₕ
        λ, (sᵢ, a, cᵣ, cₜ, cₘ, cₕ, pᵣ, pₜ, pₘ, pₕ) --> ∅

        # Non-mass action reactions
        ν_imp, ∅ => sᵢ
        ν_cat, sᵢ => nₛ * a
        νᵣ, nᵣ * a + cᵣ => pᵣ + mᵣ + pᵣ
        νₜ, nₜ * a + cₜ => pᵣ + mₜ + pₜ
        νₘ, nₘ * a + cₘ => pᵣ + mₘ + pₘ
        νₕ, nₕ * a + cₕ => pᵣ + mₕ + pₕ

        @equations begin
            λ ~ γ * (cᵣ + cₜ + cₘ + cₕ) / Mₛ
            double ~ log(2) / λ
            M ~ nᵣ * pᵣ + nₜ * pₜ + nₘ * pₘ + nₕ * pₕ + nᵣ * (cᵣ + cₜ + cₘ + cₕ)
            γ ~ γ_max * a / (Kᵧ + a)
            ν_imp ~ pₜ * νₜ_max * s / (Kₜ + s)
            ν_cat ~ pₘ * νₘ_max * sᵢ / (Kₘ + sᵢ)
            νᵣ ~ cᵣ * γ / nᵣ
            νₜ ~ cₜ * γ / nₜ
            νₘ ~ cₘ * γ / nₘ
            νₕ ~ cₕ * γ / nₕ
            wᵣ ~ wᵣ_max * a / (θᵣ + a)
            wₜ ~ wₜ_max * a / (θₜ + a)
            wₘ ~ wₘ_max * a / (θₘ + a)
            wₕ ~ wₕ_max * a / (θₕ + a) / (1 + (pₕ/Kₕ)^hₕ)
        end
    end
    return growth
end

u0 = [:mᵣ => 0.0, :mₜ => 0.0, :mₘ => 0.0, :mₕ => 0.0, :sᵢ => 0.0, :a => 1000.0, :pᵣ => 500.0, :pₜ => 0.0, :pₘ => 0.0, :pₕ => 0.0, :cᵣ => 0.0, :cₜ => 0.0, :cₘ => 0.0, :cₕ => 0.0]
u0_steady_state = [:mᵣ => 105.0, :mₜ => 16.0, :mₘ => 16.0, :mₕ => 747.0, :sᵢ => 128.0, :a => 22.0, :pᵣ => 1218.0, :pₜ => 4623.0, :pₘ => 4623.0, :pₕ => 214477.0, :cᵣ => 1049.0, :cₜ => 44.0, :cₘ => 44.0, :cₕ => 2052.0]
tspan = (0., 1500.)
ps = [:s => s,
      :γ_max => γ_max,
      :νₜ_max => νₜ_max,
      :νₘ_max => νₘ_max,
      :wᵣ_max => wᵣ_max,
      :wₕ_max => wₕ_max,
      :Kₕ => Kₕ,
      :θₜ => θₜ,
      :θₘ => θₘ,
      :θₕ => θₕ,
      :θᵣ => θᵣ,
      :k => k,
      :Mₛ => Mₛ,
      :Kᵧ => Kᵧ,
      :Kₜ => Kₜ,
      :Kₘ => Kₘ,
      :wₘ_max => wₘ_max,
      :wₜ_max => wₜ_max,
      :dₘ => dₘ,
      :hₕ => hₕ,
      :nᵣ => nᵣ,
      :nₕ => nₕ,
      :nₜ => nₜ,
      :nₘ => nₘ,
      :kᵤ => kᵤ,
      :nₛ => nₛ]

function ode_growth()
    growth = get_growth_model()
    return ODEProblem(growth, u0_steady_state, tspan, ps; structural_simplify = true)
end

function plot_ode()
    prob = ode_growth()
    sol = solve(prob, Tsit5(), maxiters=1e7, isoutofdomain=(u,p,t)->any(x->x<0,u))

    display(Plots.plot(sol, vars=[:M], xlabel="Time", ylabel="Concentration", title="Repressilator Protein Dynamics", linewidth=2))
end
