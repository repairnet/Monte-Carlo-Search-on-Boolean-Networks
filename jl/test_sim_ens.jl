#
# julia --project=BooleanNetworks.jl %
#

using BooleanNetworks

import Random
Random.seed!(1234)

ens_dir = "../data/WT_ensemble"
bnets = readdir(ens_dir, join=true)[1:30]
n = length(bnets)
println("Loading $n bnets...")
@time bns = load_bnet_ensemble(bnets)
println("Done")

