import glob
from timeit import default_timer as timer

print("~~~ init Julia..")
from juliacall import Main as jl
print("--- init Julia done")

jl.seval("using BooleanNetworks")

ens_dir = "../data/WT_ensemble"
bnets = glob.glob(f"{ens_dir}/*.bnet")[:30]

print(f"Loading {len(bnets)} bnets..")
t0 = timer()
bns = jl.load_bnet_ensemble(bnets)
print(f"Done [{timer() -t0}s]")

init_active = ["miR200", "miR203", "miR34"]
x0 = jl.zerocfg(bns[0])
for a in init_active:
    x0[bns[0].index[a]] = True

iApop = bns[0].index["Apoptosis"]
iCC = bns[0].index["CellCycleArrest"]
iInv = bns[0].index["Invasion"]
iEMT = bns[0].index["EMT"]
in_target = jl.seval(f"in_target(x) = ~x[{iApop}] & x[{iCC}] & x[{iInv}] & x[{iEMT}]")

nb_sims = 10_000
maxsteps = 300

print(f"Warmup...")
t0 = timer()
jl.fasync_ping(bns, 1, 2, x0, in_target)
print(f"Warmup done [{timer() - t0}]")

nb_sims_per = max(1, nb_sims // len(bns))

# wild type
print(f"Ping {len(bns)} BNs with {nb_sims_per} simulation each")
t0 = timer()
count = jl.fasync_ping(bns, nb_sims_per, maxsteps, x0, in_target)
print(f"Ping = {100*count / (nb_sims_per*len(bns))}% [{timer()-t0}s]")


# mutant
mutant = {"NICD": True, "p53": False}
bns_m = jl.make_mutant(bns, mutant)
print(f"Ping {len(bns)} BNs with mutant {mutant} with {nb_sims_per} simulation each")
t0 = timer()
count = jl.fasync_ping(bns_m, nb_sims_per, maxsteps, x0, in_target)
print(f"Ping = {100*count / (nb_sims_per*len(bns))}% [{timer()-t0}s]")

