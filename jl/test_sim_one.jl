#
# julia --project=BooleanNetworks.jl %
#

using BooleanNetworks

import Random
Random.seed!(1234)

filename = "../data/WT_ensemble/bn0.bnet"
bn = load_bnet(filename)

init_active = ["miR200", "miR203", "miR34"]
x0 = zerocfg(bn)
for a in init_active
    x0[bn.index[a]] = true
end

iApop = bn.index["Apoptosis"]
iCC = bn.index["CellCycleArrest"]
iInv = bn.index["Invasion"]
iEMT = bn.index["EMT"]
in_target(x) = ~x[iApop] & x[iCC] & x[iInv] & x[iEMT]

nb_sims = 10000
maxsteps = 300

# warm-up
bn_m = make_mutant(bn, zerocfg(bn), zerocfg(bn))
fasync_ping([bn,bn_m], 1, 10, x0, in_target)


println("Ping with $nb_sims simulations")
@time result = fasync_ping(bn, nb_sims, maxsteps, x0, in_target)
println("Ping = $(100*result/nb_sims)%")

mutant = Dict("NICD" => true, "p53" => false)
bn_m = make_mutant(bn, mutant)

println("Ping with $nb_sims simulations of mutant $mutant")
@time result = fasync_ping(bn_m, nb_sims, maxsteps, x0, in_target)
println("Ping = $(100*result/nb_sims)%")

println("Ping with $nb_sims * 2 simulations")
@time result = fasync_ping([bn,bn_m], nb_sims, maxsteps, x0, in_target)
println("Ping = $(100*result/(2*nb_sims))%")

