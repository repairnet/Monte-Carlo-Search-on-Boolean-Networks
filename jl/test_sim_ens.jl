#
# julia --project=BooleanNetworks.jl %
#

using BooleanNetworks
using ProgressBars

import Random
Random.seed!(1234)

ens_dir = "../data/WT_ensemble"
bnets = readdir(ens_dir, join=true)[1:30]
n = length(bnets)
println("Loading $n bnets...")
@time bns = load_bnet_ensemble(ProgressBar(bnets))
println("Done")

init_active = ["miR200", "miR203", "miR34"]
x0 = zerocfg(bns[1])
for a in init_active
    x0[bns[1].index[a]] = true
end

iApop = bns[1].index["Apoptosis"]
iCC = bns[1].index["CellCycleArrest"]
iInv = bns[1].index["Invasion"]
iEMT = bns[1].index["EMT"]
in_target(x) = ~x[iApop] & x[iCC] & x[iInv] & x[iEMT]

nb_sims = 10_000
maxsteps = 300

# warm-up (TODO: should be done at compilation..)
println("Warmup...")
@time bns_m = make_mutant(bns, Dict("Apoptosis" => true))
@time fasync_ping(bns, 1, 2, x0, in_target)
println("Warmup done.")

nb_sims_per = max(1, nb_sims ÷ length(bns))


# wild type
println("Ping $(length(bns)) BNs with $nb_sims_per simulation each")
@time count = fasync_ping(bns, nb_sims_per, maxsteps, x0, in_target)
println("Ping = $(100*count / (nb_sims_per*length(bns)))%")


# mutant
mutant = Dict("NICD" => true, "p53" => false)
bns_m = make_mutant(bns, mutant)
println("Ping with $nb_sims simulations of mutant $mutant")
@time count = fasync_ping(bns_m, nb_sims_per, maxsteps, x0, in_target)
println("Ping = $(100*count / (nb_sims_per*length(bns)))%")

