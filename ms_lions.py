import msprime
import numpy as np
import sys
import hashlib
import gzip

### FORWARD SIMULATION PARAMETERS

### Simulations carried out in msprime (v1.2.0) using a range of divergence times and rates of migration between modern and cave lions. 
### Assuming a split of an ancestral lion lineage split from an outgroup 6.19 MYA (Yuan et al 2024)
### Then splitting into modern and cave lion lineages 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75 and 2 ma. 
### A split between a North and South cave lion populations (of equal effective population size) 200 ka 
### A split between North and South modern lion lineages 70 ka
### Continuous gene-flow from the North to the South cave lion population (rate = 0.9). 
### Bi-directional gene-flow between the southern cave lion and northern modern lion populations, at rates of between 1.0e-6 to 3.0e-5, in steps of 1.0e-6 (total of 30 steps)
### Modern/cave gene-flow modelled as “pulses”, occurring 100 years before each sampled cave lion genome, lasting for 1 generation
### Generation time of five years 
### Estimates of effective population size were taken from the above PSMC analysis.

# ---------------------- Parse Arguments ---------------------- #
if len(sys.argv) < 6:
    print("Usage: python script.py <chromosome_name> <chromosome_length> <uniform_mig_rate> <pulse_offset_years> <cave_lion_migration_rate>")
    sys.exit(1)

chrom_name = sys.argv[1]
sequence_length = float(sys.argv[2])  # chr length
uniform_mig_rate = float(sys.argv[3])  # modern <-> cave migration rate
pulse_offset_years = float(sys.argv[4])  # how long before the sample the pulse happens
cave_lion_migration_rate = float(sys.argv[5])  # migration rate PopCaveN -> PopCaveS

## Run as e.g: python ms_Lions.py Mb50_5e-6_100 50000000 0.000005 100 0.900

# ---------------------- Deterministic Random Seed ---------------------- #
def make_seed(name):
    return int(hashlib.md5(name.encode()).hexdigest(), 16) % (2**32)

seed = make_seed(chrom_name)
print(f"Using random seed {seed} for chromosome {chrom_name}")
print(f"Migration pulses: rate={uniform_mig_rate}, 1 generation long, "
      f"occurring {pulse_offset_years} years before each ancient genome is sampled.")
print("Migration rate units: per-generation probability that a lineage in one population migrates to the other population.")
print(f"Cave lion continuous migration rate (PopCaveN → PopCaveS): {cave_lion_migration_rate}")
print("Modern lion continuous migration rate (PopModN ↔ PopModS): 0.000 (NO MIGRATION)")

# ---------------------- Parameters ---------------------- #
Ne_Cave = 31558
Ne_Mod = 27029
Ne_anc = 86679
Ne_out = 20000

mutation_rate = 4.5e-9
recomb_rate = 1e-8
sample_size_Mod = 21
outgroup_split_time = 1238000  # generations
generation_time = 5  # years

# Divergence times
divergence_times = np.arange(50000, 450000, 50000)  # in generations

# ---------------------- Ancient Samples (PopCave) ---------------------- #
popCave_metadata = {
    "Austria17KYApseu": 17000,
    "Siberia44KYApseu": 44000,
    "Siberia30KYApseu": 30000,
    "SiberiaA33KYApseu": 33000,
    "NSI52KYApseu": 52000,
    "Siberia47KYApseu": 47000,
    "Siberia22KYApseu": 22000,
    "Siberia148KYApseu": 148000,
    "Siberia64KYApseu": 64000,
    "CEAsia20KYA": 20000,
}

# ---------------------- PopMod Samples ---------------------- #
popModN_names = [
    "1916pseu", "BarbaryRB673dip", "IranRB717pseu", "BarbaryRB687pseu",
    "BarbaryRB713pseu", "BarbaryNHMRB753", "SeneRB691pseu",
    "SeneRB710pseu", "GabonRB709pseu", "India1pseu", "India2pseu"
]

popModS_names = [
    "Botswana1pseu", "Botswana2pseu", "CapeRB668pseu", "CapeRB720pseu",
    "Pantheraleo1dip", "Pantheraleo2dip", "RSARB703pseu", "SudanRB742pseu",
    "Tanzania1pseu", "Tanzania2pseu"
]

# ---------------------- Main Simulation Loop ---------------------- #
for base_div_time in divergence_times:
    print(f"\nSimulating for divergence time {base_div_time} generations "
          f"({base_div_time * generation_time:,} years)")

    demography = msprime.Demography()

    # Define populations
    demography.add_population(name="PopCaveAnc", initial_size=Ne_Cave)

    # At split, PopCaveN = 90% of Ne_Cave, PopCaveS = 10% of Ne_Cave
    demography.add_population(name="PopCaveS", initial_size=int(0.1 * Ne_Cave))
    demography.add_population(name="PopCaveN", initial_size=int(0.9 * Ne_Cave))

    demography.add_population(name="PopModAnc", initial_size=Ne_Mod)
    demography.add_population(name="PopModN", initial_size=0.5 * Ne_Mod)
    demography.add_population(name="PopModS", initial_size=0.5 * Ne_Mod)

    demography.add_population(name="anc", initial_size=Ne_anc)
    demography.add_population(name="outgroup", initial_size=Ne_out)

    print(f"At split: PopCaveN size = {int(0.9 * Ne_Cave)}, PopCaveS size = {int(0.1 * Ne_Cave)}")
    print(f"At split: PopModN size = {int(0.5 * Ne_Mod)}, PopModS size = {int(0.5 * Ne_Mod)}")

    # Continuous migration
    demography.set_migration_rate(source="PopCaveN", dest="PopCaveS", rate=cave_lion_migration_rate)
    demography.set_symmetric_migration_rate(["PopModS","PopModN"], rate=0.000)

    # ---------------------- 1-gen migration pulses ---------------------- #
    pulse_offset_gens = int(pulse_offset_years // generation_time)

    for target_years in [64000, 52000, 47000, 44000, 33000, 30000, 22000, 20000, 17000]:
        sample_time = target_years // generation_time

        pulse_end = sample_time + pulse_offset_gens
        pulse_start = pulse_end + 1   # exactly 1 generation

        print(f"Migration pulse before {target_years} years ago: "
              f"{pulse_start}-{pulse_end} generations ago "
              f"(1 generation long, ending {pulse_end*generation_time} years before sampling) "
              f"at rate {uniform_mig_rate}")

        demography.add_symmetric_migration_rate_change(
            time=pulse_start,
            populations=["PopCaveS","PopModN"],
            rate=uniform_mig_rate
        )
        demography.add_symmetric_migration_rate_change(
            time=pulse_end,
            populations=["PopCaveS","PopModN"],
            rate=0
        )

    # ---------------------- Splits ---------------------- #
    demography.add_population_split(time=200000//generation_time, derived=["PopCaveS","PopCaveN"], ancestral="PopCaveAnc")
    demography.add_population_split(time=70000//generation_time, derived=["PopModN","PopModS"], ancestral="PopModAnc")
    demography.add_population_split(time=base_div_time, derived=["PopCaveAnc","PopModAnc"], ancestral="anc")
    demography.add_population_split(time=outgroup_split_time, derived=["anc"], ancestral="outgroup")

    demography.sort_events()

    # ---------------------- Samples ---------------------- #
    samples = [
        msprime.SampleSet(1, population="PopCaveS", ploidy=2, time=age_years//generation_time)
        if sample_name not in ["Siberia148KYApseu","NSI52KYApseu"] else
        msprime.SampleSet(1, population="PopCaveN", ploidy=2, time=age_years//generation_time)
        for sample_name, age_years in popCave_metadata.items()
    ]

    samples += [msprime.SampleSet(1, population="PopModN", ploidy=2) for _ in popModN_names]
    samples += [msprime.SampleSet(1, population="PopModS", ploidy=2) for _ in popModS_names]
    samples.append(msprime.SampleSet(1, population="outgroup", ploidy=2))

    # ---------------------- Run Simulation ---------------------- #
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recomb_rate,
        random_seed=seed
    )
    mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seed)

    # ---------------------- Thinning (optional) ---------------------- #
    thinning_factor = 1
    sites = list(mts.sites())
    keep_every = int(np.ceil(thinning_factor))
    keep_indices = np.arange(0, len(sites), keep_every)
    keep_sites = [sites[i].id for i in keep_indices]
    mts_thinned = mts.delete_sites([s.id for s in sites if s.id not in keep_sites])

    # ---------------------- Write VCF ---------------------- #
    sample_names = list(popCave_metadata.keys()) + popModN_names + popModS_names + ["outgroup"]
    vcf_file = f"{chrom_name}_ALL_ancients_div_{base_div_time}_gens.vcf.gz"
    with gzip.open(vcf_file, "wt") as vcf_output:
        mts_thinned.write_vcf(vcf_output, individual_names=sample_names)

    print(f"Wrote {vcf_file} with {len(sample_names)} samples.")
