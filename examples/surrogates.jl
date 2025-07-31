using Plots
default(c=:matter, legend=false, xlabel="x", ylabel="y")
using Surrogates
using Flux
using SurrogatesFlux

n_samples = 1000

lb = 1.0
ub = 6.0
objective = z -> z + 1
xys = sample(n_samples, lb, ub, SobolSample())
zs = objective.(xys)
plt = scatter(xys, zs, label="Samples", color=:blue)
if length(lb) == 1
    xys = [[x] for x in xys]
end
model = Chain(
    Dense(length(lb), 32, relu),
    Dense(32, 32, relu),
    Dense(32, length(lb)),
    first
)
neural = NeuralSurrogate(xys, zs, lb, ub, model=model, n_epochs=1000)
plot!(plt, lb:0.01:ub, [neural(x) for x in lb:0.01:ub], label="Surrogate (untrained)", color=:red)
for (x, z) in zip(xys, zs)
    update!(neural, x, z)
end
plot!(plt, lb:0.01:ub, [neural(x) for x in lb:0.01:ub], label="Surrogate (trained)", color=:green)
