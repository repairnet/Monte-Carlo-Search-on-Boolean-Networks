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

nb_sims = 10000
maxsteps = 300
println("Performing $nb_sims simulations of at most $maxsteps from fixed initial configuration")
@time result = fasync_simulations(bn, outputs, nb_sims, maxsteps, x0)

ratios = outputs_ratios(result, outputs, bn)
println(ratios)

mutated = zerocfg(bn)
mutant = zerocfg(bn)

fasync_simulations(bn, mutated, mutant, outputs, 1, 10, x0)

println("Performing $nb_sims simulations with no mutant of at most $maxsteps from fixed initial configuration")
@time result = fasync_simulations(bn, mutated, mutant, outputs, nb_sims, maxsteps, x0)
ratios = outputs_ratios(result, outputs, bn)
println(ratios)

mutated[bn.index["NICD"]] = true
mutated[bn.index["p53"]] = true
mutant[bn.index["NICD"]] = true
mutant[bn.index["p53"]] = false
println("Performing $nb_sims simulations with mutant $(bn.nodes[mutated]) of at most $maxsteps from fixed initial configuration")
@time result = fasync_simulations(bn, mutated, mutant, outputs, nb_sims, maxsteps, x0)
ratios = outputs_ratios(result, outputs, bn)
println(ratios)

using Profile
@profile fasync_simulations(bn, mutated, mutant, outputs, nb_sims, maxsteps, x0)
Profile.print()

