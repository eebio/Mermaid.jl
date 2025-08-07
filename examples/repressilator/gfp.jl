using Catalyst
using DifferentialEquations, StochasticDiffEq
using Plots

function get_repressilator()
    # Constants
    eff = 20
    n = 2
    KM = 40
    tau_mRNA = 2
    tau_prot = 10
    ps_a = 0.5
    ps_0 = 5e-4

    # Derived
    t_ave = tau_mRNA / log(2)
    beta = tau_mRNA / tau_prot
    k_tl = eff / t_ave
    a_tr = (ps_a - ps_0) * 60
    a0_tr = ps_0 * 60
    kd_prot = log(2) / tau_prot
    kd_mRNA = log(2) / tau_mRNA
    alpha = a_tr * eff * tau_prot / (KM * log(2))

    repressilator = @reaction_network repressilator begin
        hillr(P₃, a_tr, KM, n) + a0_tr, ∅ --> m₁
        hillr(P₁, a_tr, KM, n) + a0_tr, ∅ --> m₂
        hillr(P₂, a_tr, KM, n) + a0_tr, ∅ --> m₃
        kd_mRNA, m₁ --> ∅
        kd_mRNA, m₂ --> ∅
        kd_mRNA, m₃ --> ∅
        k_tl, m₁ --> m₁ + P₁
        k_tl, m₂ --> m₂ + P₂
        k_tl, m₃ --> m₃ + P₃
        kd_prot, P₁ --> ∅
        kd_prot, P₂ --> ∅
        kd_prot, P₃ --> ∅
        hillr(P₁, a_tr, KM, n) + a0_tr, ∅ --> m_gfp
        kd_mRNA, m_gfp --> ∅
        k_tl, m_gfp --> m_gfp + gfp
        kd_prot, gfp --> ∅
    end
    return repressilator
end

u0 = [:P₁ => 150.0, :P₂ => 140.0, :P₃ => 130.0, :m₁ => 10.0, :m₂ => 12.0, :m₃ => 11.0, :m_gfp => 0.0, :gfp => 0.0]
tspan = (0., 1500.)
ps = [:k_tl => k_tl, :KM => KM, :a0_tr => a0_tr, :a_tr => a_tr, :kd_prot => kd_prot, :n => n, :kd_mRNA => kd_mRNA]

function ode_repressilator()
    repressilator = get_repressilator()
    return ODEProblem(repressilator, u0, tspan, ps)
end

function plot_ode()
    prob = ode_repressilator()
    sol = solve(prob, Tsit5())

    plot(sol, vars=[:gfp], xlabel="Time", ylabel="Concentration", label="gfp - Cell 1", title="Repressilator Protein Dynamics", linewidth=2)
    sol = solve(prob, Tsit5())
    display(plot!(sol, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))
end

function sde_repressilator()
    repressilator = get_repressilator()
    return SDEProblem(repressilator, u0, tspan, ps)
end

function plot_sde()
    sde = sde_repressilator()
    sol_sde = solve(sde, EM(), dt=0.1)

    plot(sol_sde, vars=[:gfp], xlabel="Time", ylabel="Concentration", label="gfp - Cell 1", title="Repressilator Protein Dynamics (Stochastic)", linewidth=2)
    sol_sde = solve(sde, EM(), dt=0.1)
    display(plot!(sol_sde, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))
end
