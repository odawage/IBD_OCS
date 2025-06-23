
import numpy as np
import math
import msprime
import gzip
import stdpopsim
import tskit 
import os
import sys


path = sys.argv[1]
Chr = sys.argv[2]
sample_set = sys.argv[3] # number of individuals to create
NE = sys.argv[4] #NE


full_path = os.path.join(path, Chr + '.vcf.gz')
print(full_path)


# Inputs 
# NE = 1000
#sample_set = 10
ploidy = 2
mutation_rate=float(2.2e-9) # 2.2 x 10(-9) 0.94e-8
recombination_rate = 1e-8

# genome length and chromosome length 
genome_length = 3e9 # 3Gbp pig genome is 2.5 - 2.7Gbp long, or 2 700 000 000
NChr = 30
Chr = int(Chr)
# Since it's a non real animals we do equal length chromosomes 
chr_length = genome_length/NChr 

# Then we set up the model. This should give us a constant NE pop. 
# should see if the coalecens matches up 
model = msprime.Demography.isolated_model(
    initial_size=[NE], 
    growth_rate=[0]
)
# Constant recomb. across the chromosome 
recombination_map = msprime.RateMap.uniform(
    sequence_length=int(chr_length), 
    rate=recombination_rate
)
# generate seeds to sample from
rng = np.random.RandomState(42)
seeds = rng.randint(1, 2**31, size=(NChr, 2))
ancestry_seed, mutation_seed = seeds[NChr-1]

# then we simulate the ancestry 
ts = msprime.sim_ancestry(
    samples=int(sample_set),
    ploidy=ploidy,
    recombination_rate=recombination_map,
    random_seed=int(ancestry_seed),
    demography=model,
    end_time=None
)

# for future self: 
#   print(ts.first().draw_text()) A way to draw the tree
#   population_size=int(NE),

# then add mutations to the simulation 
ts = msprime.sim_mutations(
    ts,
    end_time=None,
    random_seed=mutation_seed,
    rate=mutation_rate
)

# Write the results 
with gzip.open(full_path, "wt") as f:
    ts.write_vcf(f,contig_id=Chr, position_transform = lambda x: np.fmax(1, x))

