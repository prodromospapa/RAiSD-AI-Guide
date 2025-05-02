#!/usr/bin/env python
import stdpopsim
from matplotlib import pyplot as plt
import argparse
import os
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# Command line arguments
parser = argparse.ArgumentParser(description="Simulate population genetics data.")
parser.add_argument("-dir", type=str, help="Output directory")
parser.add_argument("-ntrain", type=int, help="Number of training samples")
parser.add_argument("-ntest", type=int, help="Number of test samples")
parser.add_argument("-train_pos", type=float, help="train sweep position")
parser.add_argument("-test_pos", type=float, help="test sweep position")
parser.add_argument("-species", type=str, help="Species to simulate")
parser.add_argument("-chr", type=str, help="Chromosome to simulate")
parser.add_argument("-max_samples", type=int, help="Max samples per population")
args = parser.parse_args()

# Create output directories
dir=args.dir
if not dir.endswith("/"):
    dir += "/"
species_id = args.species
ntrain = args.ntrain
ntest = args.ntest
train_pos = args.train_pos
test_pos = args.test_pos
max_samples_per_pop = args.max_samples
chr = args.chr

species = stdpopsim.get_species(species_id)
Ne = species.population_size

contig = species.get_contig(chr)
mu = contig.mutation_rate
recombination_rate = contig.recombination_map.mean_rate
ploidy=contig.ploidy
L = int(contig.length)
engine = stdpopsim.get_engine("slim")

model = stdpopsim.PiecewiseConstantSize(Ne)

def selective_sweep_model(pos):
    id = f"hard_sweep_{np.random.randint(0, 100000)}"
    contig.add_single_site(
        id=id,
        coordinate=pos)
    
    extended_events = stdpopsim.selective_sweep(
        single_site_id=id,
        population="pop_0",
        selection_coeff=0.5,
        mutation_generation_ago=1000,
    )
    return extended_events


def simulate_replicates(samples,pos=False):
    if pos:
        extended_events = selective_sweep_model(pos)
    else:
        extended_events = None
    return engine.simulate(
    model,
    contig,
    {"pop_0": samples},
    extended_events=extended_events,
    slim_scaling_factor=10,
    slim_burn_in=0.1
) 

# Function to save multiple simulations into one file
def save2vcf(samples,filepath,pos=False,max_samples_per_pop=max_samples_per_pop):
    i=0
    file_list = []
    while samples >0:
        i+=1
        part = min(samples,max_samples_per_pop)
        file_list.append(filepath+str(i))
        with open(filepath+str(i), "w") as f: 
            ts = simulate_replicates(part,pos)
            ts.write_vcf(f,contig_id=i)
        samples -= part
    # merge files
    os.system(f"bcftools concat -o {filepath} " + " ".join(file_list) + " >/dev/null 2>&1")
    # remove temp files
    for file in file_list:
        os.remove(file)

os.makedirs(dir, exist_ok=True)
# Simulate neutral and hard sweep and save to vcf
save2vcf(ntest, f"{dir}neutral.vcf")
print("Neutral simulation done")
save2vcf(ntrain, f"{dir}train_sweep.vcf",train_pos)
print("Sweep train simulation done")
save2vcf(ntest,f"{dir}test_sweep.vcf",test_pos)
print("Sweep test simulation done")