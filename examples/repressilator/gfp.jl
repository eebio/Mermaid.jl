module Repressilator

    using OrdinaryDiffEq
    using StochasticDiffEq
    using Catalyst
    using JumpProcesses
    using Plots
    using Random

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

    repressilator = @reaction_network begin
        @species gr(t) V(t)
        a*gr*V*(hillr(P₃/V, a_tr, KM, n) + a0_tr), ∅ --> m₁
        a*gr*V*(hillr(P₁/V, a_tr, KM, n) + a0_tr), ∅ --> m₂
        a*gr*V*(hillr(P₂/V, a_tr, KM, n) + a0_tr), ∅ --> m₃
        b*gr*kd_mRNA, m₁ --> ∅
        b*gr*kd_mRNA, m₂ --> ∅
        b*gr*kd_mRNA, m₃ --> ∅
        a*gr*k_tl, m₁ --> m₁ + P₁
        a*gr*k_tl, m₂ --> m₂ + P₂
        a*gr*k_tl, m₃ --> m₃ + P₃
        b*gr*kd_prot, P₁ --> ∅
        b*gr*kd_prot, P₂ --> ∅
        b*gr*kd_prot, P₃ --> ∅
        a*gr*V*(hillr(P₁/V, a_tr, KM, n) + a0_tr), ∅ --> m_gfp
        b*gr*kd_mRNA, m_gfp --> ∅
        b*gr*k_tl, m_gfp --> m_gfp + gfp
        b*gr*kd_prot, gfp --> ∅
    end

    u0 = [:P₁ => 10.0, :P₂ => 10.0, :P₃ => 500.0, :m₁ => 50.0, :m₂ => 5.0, :m₃ => 5.0, :m_gfp => 0.0, :gfp => 0.0, :gr => 1.0, :V => 1.0]
    tspan = (0., 750.)
    ps = [:k_tl => k_tl, :KM => KM, :a0_tr => a0_tr, :a_tr => a_tr, :kd_prot => kd_prot, :n => n, :kd_mRNA => kd_mRNA, :a => 10.0, :b => 40.0]

    function ode_repressilator()
        return ODEProblem(repressilator, u0, tspan, ps)
    end

    function plot_ode()
        prob = ode_repressilator()
        sol = solve(prob, Tsit5())

        Plots.plot(sol, vars=[:gfp], xlabel="Time", ylabel="Concentration", label="gfp - Cell 1", title="Repressilator Protein Dynamics", linewidth=2)
        sol = solve(prob, Tsit5())
        display(Plots.plot!(sol, vars=[:gfp], label="gfp - Cell 2", linestyle=:dash, linewidth=2))
    end

    function jump_repressilator()
    u0_integers = [:P₁ => 10, :P₂ => 10, :P₃ => 500, :m₁ => 50, :m₂ => 5,
        :m₃ => 5, :m_gfp => 0, :gfp => 0, :gr => 1.0, :V => 1.0]

        prob = DiscreteProblem(repressilator, u0_integers, tspan, ps)
        jprob = JumpProblem(repressilator, prob, Direct(); save_positions = (false, false))
        return jprob
    end

    function plot_jump()
        prob = jump_repressilator()
        sol = solve(prob)

        Plots.plot(sol, vars = [:gfp], xlabel = "Time", ylabel = "Concentration",
            label = "gfp - Cell 1", title = "Repressilator Protein Dynamics", linewidth = 2)
        sol = solve(prob)
        display(Plots.plot!(
            sol, vars = [:gfp], label = "gfp - Cell 2", linestyle = :dash, linewidth = 2))
    end

    function sde_repressilator()
        return SDEProblem(repressilator, u0, tspan, ps)
    end

    function plot_sde()
        sde = sde_repressilator()
        sol_sde = solve(sde, EM(), dt = 0.1)

        Plots.plot(sol_sde, vars = [:gfp], xlabel = "Time",
            ylabel = "Concentration", label = "gfp - Cell 1",
            title = "Repressilator Protein Dynamics (Stochastic)", linewidth = 2)
        sol_sde = solve(sde, EM(), dt = 0.1)
        display(Plots.plot!(
            sol_sde, vars = [:gfp], label = "gfp - Cell 2", linestyle = :dash, linewidth = 2))
    end
end
using .Repressilator
