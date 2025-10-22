#!/bin/bash -l
#SBATCH --nodes=1 
#SBATCH  --nodelist=cbsubscb18            
#SBATCH --ntasks=1               
#SBATCH --mem=10000               
#SBATCH --time=0-24:00:00       
#SBATCH --partition=regular         
#SBATCH --chdir=/NFS4/storage/Jinmin/slurm
#SBATCH --job-name=convert        
#SBATCH --output=convert.out.%j  
#SBATCH --mail-user=jl4566@cornell.edu
#SBATCH --mail-type=ALL   
#SBATCH --array=0-9

# ---- PARAMETERS ----

source /fs/cbsubscb18/storage/Jinmin/myenv/bin/activate

name="Ne2e+04_r1e-08_100_1e+06"
singer_bin="/NFS4/storage/Jinmin/SINGER/releases/singer_master"
convert_bin="/fs/cbsubscb18/storage/Jinmin/singer2/SINGER/releases/convert_to_tskit"

# Each array job handles 10 replicates
start_idx=$(( SLURM_ARRAY_TASK_ID * 10 ))
end_idx=$(( start_idx + 9 ))

echo "Node $(hostname) handling replicates ${start_idx}-${end_idx}"

for i in $(seq ${start_idx} ${end_idx}); do
    echo "Processing replicate $i ..."

    vcf_in="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}.vcf.gz"
    vcf_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}.vcf"
    singer_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/singer_temp/${name}/${i}/${name}_${i}"
    tskit_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/tree/${name}/${i}/singer_inferred/${name}"

    mkdir -p "$(dirname "${vcf_out}")" "$(dirname "${singer_out}")" "$(dirname "${tskit_out}")"

    # # ---- Convert gzipped VCF ----
    # echo "Converting ${vcf_in} -> ${vcf_out}"
    # gunzip -c "${vcf_in}" > "${vcf_out}"

    # # ---- Run SINGER ----
    # echo "Running SINGER for ${vcf_out}"
    # "${singer_bin}" \
    #     -Ne 2e4 \
    #     -m 1.25e-8 \
    #     -vcf "/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}" \
    #     -output "${singer_out}" \
    #     -start 0 \
    #     -end 1e6 \
    #     -thin 50 \
    #     -n 50 \
    #     -polar 0.99

    # ---- Convert to TSKIT ----
    echo "Converting SINGER output to tskit format"
    "${convert_bin}" \
        -input "${singer_out}" \
        -output "${tskit_out}" \
        -start 30 \
        -end 50 \
        -step 1

    echo "Finished replicate $i"
    echo "---------------------------------------------"
done

echo "Node $(hostname) finished replicates ${start_idx}-${end_idx}."
