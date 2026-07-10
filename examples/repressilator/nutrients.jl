module Nutrients
    using OrdinaryDiffEq
    using Catalyst
    function get_nutrient_prob()
        rs = @reaction_network begin
            @species nut_rate(t)
            100*nut_rate, nut => ∅ # Dont want this to be mass action
        end
        diffusion = @transport_reaction D nut
        space = CartesianGrid((100,100))
        dsrs = DiscreteSpaceReactionSystem(rs, [diffusion], space)
        u0 = [:nut_rate => 0.0, :nut => 1.0]
        tspan = (0.0, 15.0)
        ps = [:D => 1e-1]
        prob = ODEProblem(dsrs, u0, tspan, ps)
        return prob, dsrs
    end
end
