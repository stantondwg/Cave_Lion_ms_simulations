# Cave Lion msprime simulations

This repository contains a single Python script, `ms_lions.py`, which uses **msprime** to simulate genetic variation under a demographic model involving **modern lions**, **cave lions**, and an **outgroup**.

## What the script does (high level)

For a range of divergence times between modern and cave lion lineages, the script:

1. Builds an `msprime.Demography()` model with multiple populations (cave lion North/South, modern lion North/South, ancestral populations, and an outgroup).
2. Adds population split events at specified times (in generations).
3. Applies:
   - **Continuous migration** from North cave lions → South cave lions.
   - **Short (1-generation) migration pulses** between South cave lions and North modern lions that occur a fixed number of years before each ancient cave lion sample time.
4. Simulates ancestry with recombination (`msprime.sim_ancestry`) and then overlays mutations (`msprime.sim_mutations`).
5. Writes a gzipped VCF for each divergence-time replicate.

## Script inputs

The script is run from the command line and expects **five** arguments:

```bash
python ms_lions.py <chromosome_name> <chromosome_length> <uniform_mig_rate> <pulse_offset_years> <cave_lion_migration_rate>
```

- `chromosome_name`: A label used for output naming and to derive a deterministic random seed.
- `chromosome_length`: Sequence length to simulate (e.g., `50000000` for 50 Mb).
- `uniform_mig_rate`: Migration rate used during the 1-generation pulses between **PopCaveS** and **PopModN**.
- `pulse_offset_years`: How many years before each ancient genome’s sampling time the pulse occurs.
- `cave_lion_migration_rate`: Continuous migration rate from **PopCaveN → PopCaveS**.

Example:

```bash
python ms_lions.py Mb50_5e-6_100 50000000 0.000005 100 0.900
```

## Outputs

For each divergence time tested, the script writes a gzipped VCF named like:

```
<chromosome_name>_ALL_ancients_div_<divergence_time>_gens.vcf.gz
```

These VCFs contain genotypes for:

- Multiple **ancient cave lion** individuals sampled at different times (years before present)
- **Modern lion** samples split into North and South groups
- A single **outgroup** sample
