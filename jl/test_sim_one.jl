#
# julia --project=BooleanNetworks.jl %
#

using BooleanNetworks

import Random
Random.seed!(1234)

filename = "../data/WT_ensemble/bn0.bnet"
@time bn = load_bnet(filename)
init_active = ["miR200", "miR203", "miR34"]
outputs = ["Apoptosis"; "CellCycleArrest"; "Invasion"; "EMT"]
outputs = [bn.index[node] for node in outputs]

x0 = zerocfg(bn)
for a in init_active
    x0[bn.index[a]] = true
end

# warm-up
fasync_simulations(bn, outputs, 1, 10, x0)

nb_sims = 100000
maxsteps = 500
println("Performing $nb_sims simulations of at most $maxsteps from fixed initial configuration")
@time result = fasync_simulations(bn, outputs, nb_sims, maxsteps, x0)

ratios = outputs_ratios(result, outputs, bn)
println(ratios)
