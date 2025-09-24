import msprime
import tskit
import math

# Times are provided in years, so we convert into generations.
generation_time = 25
T_OOA = 21.2e3 / generation_time
T_AMH = 140e3 / generation_time
T_ANC = 220e3 / generation_time
# We need to work out the starting population sizes based on
# the growth rates provided for these two populations
r_CEU = 0.004
r_CHB = 0.0055
N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
N_CHB = 510 / math.exp(-r_CHB * T_OOA)

demography = msprime.Demography()
# This is the "trunk" population that we merge other populations into
demography.add_population(
    name="YRI",
    description="Africa",
    initial_size=12300,
    # NOTE: we have to set this flag if we have a population that's
    # ancestral in a population split *and* is contemporary population.
    initially_active=True,
)
demography.add_population(
    name="CEU",
    description="European",
    initial_size=N_CEU,
    growth_rate=r_CEU,
)
demography.add_population(
    name="CHB",
    description="East Asian",
    initial_size=N_CHB,
    growth_rate=r_CHB,
)
demography.add_population(
    name="OOA",
    description="Bottleneck out-of-Africa population",
    initial_size=2100,
)

# Set the migration rates between extant populations
demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

demography.add_population_split(
    time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
)
demography.add_symmetric_migration_rate_change(
    time=T_OOA, populations=["YRI", "OOA"], rate=25e-5
)
demography.add_population_split(
    time=T_AMH, derived=["OOA"], ancestral="YRI"
)
demography.add_population_parameters_change(
    time=T_ANC, population="YRI", initial_size=7300
)
#demography.debug()


import msprime
import tskit
import math
def simple_OOA_sim(n,lens=1e6):
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_OOA = 21.2e3 / generation_time
    T_AMH = 140e3 / generation_time
    T_ANC = 220e3 / generation_time
    # We need to work out the starting population sizes based on
    # the growth rates provided for these two populations
    r_CEU = 0.004
    r_CHB = 0.0055
    N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
    N_CHB = 510 / math.exp(-r_CHB * T_OOA)

    demography = msprime.Demography()
    # This is the "trunk" population that we merge other populations into
    demography.add_population(
        name="YRI",
        description="Africa",
        initial_size=12300,
        # NOTE: we have to set this flag if we have a population that's
        # ancestral in a population split *and* is contemporary population.
        initially_active=True,
    )
    demography.add_population(
        name="CEU",
        description="European",
        initial_size=N_CEU,
        growth_rate=r_CEU,
    )
    demography.add_population(
        name="CHB",
        description="East Asian",
        initial_size=N_CHB,
        growth_rate=r_CHB,
    )
    demography.add_population(
        name="OOA",
        description="Bottleneck out-of-Africa population",
        initial_size=2100,
    )

    # Set the migration rates between extant populations
    demography.set_symmetric_migration_rate(["CEU", "CHB"], 9.6e-5)
    demography.set_symmetric_migration_rate(["YRI", "CHB"], 1.9e-5)
    demography.set_symmetric_migration_rate(["YRI", "CEU"], 3e-5)

    demography.add_population_split(
        time=T_OOA, derived=["CEU", "CHB"], ancestral="OOA"
    )
    demography.add_symmetric_migration_rate_change(
        time=T_OOA, populations=["YRI", "OOA"], rate=25e-5
    )
    demography.add_population_split(
        time=T_AMH, derived=["OOA"], ancestral="YRI"
    )
    demography.add_population_parameters_change(
        time=T_ANC, population="YRI", initial_size=7300
        )
    #demography.debug()
    new_ts = msprime.sim_ancestry(
        samples={"YRI":n,"CEU":n,"CHB":n},
        recombination_rate=1e-8,
        sequence_length=lens,
        discrete_genome=False,
        demography=demography,
        record_full_arg=False,
        #record_migrations=True,
        #model="smc_prime",
        #additional_nodes=additional_nodes,
        #coalescing_segments_only=False,
    )
    new_ts = msprime.sim_mutations(new_ts, rate=1e-8)
    return new_ts