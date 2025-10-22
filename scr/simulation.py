import stdpopsim
import os
import tskit
import msprime
import numpy as np
import tsinfer
import zarr
import demesdraw
import subprocess
import pyfaidx
import tsdate
import matplotlib.pyplot as plt
from IPython.display import display, SVG
from Bio import bgzf
#display(SVG(ts.draw_svg()))
from OOA_sim import simple_OOA_sim

def dump_ts(ts,save_dir="./simulations/sim1/",vcz_dir="./simulations/sim1/",name="test",contig_id=0):
    os.makedirs(save_dir, exist_ok=True)
    ts.dump(save_dir+name+".trees")
    np.save(f"{vcz_dir+name}-AA.npy", [s.ancestral_state for s in ts.sites()])
    vcf_name = f"{vcz_dir+name}.vcf.gz"
    with bgzf.open(vcf_name, "wt") as f:
        ts.write_vcf(f,contig_id=contig_id)
    subprocess.run(["tabix", vcf_name])
    ret = subprocess.run(
        ["python", "-m", "bio2zarr", "vcf2zarr", "convert", "--force", vcf_name, f"{vcz_dir+name}.vcz"],
        stderr = None,
    )
    if ret.returncode != 0: 
        print(f"[ERROR] bio2zarr failed for {name} with code {ret.returncode}") 
    else: 
        print(f"[OK] bio2zarr finished for {name}")
def infer_ts(vcz_full_name,tree_full_name,rate=1e-8):
    #load zarr file, ancestral states and infer
    ancestral_states = np.load(f"{vcz_full_name}-AA.npy")
    vdata = tsinfer.VariantData(f"{vcz_full_name}.vcz", ancestral_states)
    inferred_ts=tsinfer.infer(vdata,recombination_rate=1e-8)
    inferred_ts.dump(f"{tree_full_name}-inferred.trees")
    print("tsinfer done")
    print("Running tsdate")
    simplified_ts = tsdate.preprocess_ts(inferred_ts)

    redated_ts = tsdate.date(simplified_ts, mutation_rate=rate)

    out_file = f"{tree_full_name}-inferred.tsdate.trees"
    redated_ts.dump(out_file)


def simulate_dump(model,number_individuals,num_replicates,tree_dir,vcz_dir,sequence_length=1e6):
    if isinstance(model, float):
        name = f"Ne{model:.0e}_{2*number_individuals}_{sequence_length:.0e}"
        
        tree_base_path = os.path.join(tree_dir, name)
        vcz_base_path=os.path.join(vcz_dir, name)
        os.makedirs(tree_base_path, exist_ok=True)
        os.makedirs(vcz_base_path, exist_ok=True)
        for i in range(num_replicates):
            tree_replicate_dir = os.path.join(tree_base_path, str(i))+"/"
            vcz_replicate_dir = os.path.join(vcz_base_path, str(i))+"/"
            os.makedirs(tree_replicate_dir, exist_ok=True)
            os.makedirs(vcz_replicate_dir, exist_ok=True)
            ts=msprime.sim_ancestry(samples=number_individuals,population_size=model,sequence_length=sequence_length,recombination_rate=1e-8)
            ts=msprime.sim_mutations(ts,rate=1e-8)
            dump_ts(ts,save_dir=tree_replicate_dir,vcz_dir=vcz_replicate_dir,name=name+"_"+str(i))
            infer_ts(vcz_full_name=vcz_replicate_dir+name+"_"+str(i),tree_full_name=tree_replicate_dir+name+"_"+str(i))
    else:
        assert(isinstance(model, str))
        assert(model=="ooa")
        name=f"ooa_{6*number_individuals}_{sequence_length:.0e}"
    
        tree_base_path = os.path.join(tree_dir, name)
        vcz_base_path=os.path.join(vcz_dir, name)
        os.makedirs(tree_base_path, exist_ok=True)
        os.makedirs(vcz_base_path, exist_ok=True)
        for i in range(num_replicates):
            tree_replicate_dir = os.path.join(tree_base_path, str(i))+"/"
            vcz_replicate_dir = os.path.join(vcz_base_path, str(i))+"/"
            os.makedirs(tree_replicate_dir, exist_ok=True)
            os.makedirs(vcz_replicate_dir, exist_ok=True)
            ts=simple_OOA_sim(number_individuals,sequence_length)
            dump_ts(ts,save_dir=tree_replicate_dir,vcz_dir=vcz_replicate_dir,name=name+"_"+str(i))
            infer_ts(vcz_full_name=vcz_replicate_dir+name+"_"+str(i),tree_full_name=tree_replicate_dir+name+"_"+str(i))



