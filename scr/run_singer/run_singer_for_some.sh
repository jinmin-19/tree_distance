#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --nodelist=cbsubscb18   
#SBATCH --ntasks=1
#SBATCH --mem=10000
#SBATCH --time=0-24:00:00
#SBATCH --partition=regular
#SBATCH --chdir=/NFS4/storage/Jinmin/slurm
#SBATCH --job-name=runsinger_missing_all
#SBATCH --output=runsinger_missing_all.out.%A_%a
#SBATCH --mail-user=jl4566@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --array=0-19  

source /fs/cbsubscb18/storage/Jinmin/myenv/bin/activate

singer_bin="/NFS4/storage/Jinmin/SINGER/releases/singer_master"
convert_bin="/fs/cbsubscb18/storage/Jinmin/singer2/SINGER/releases/convert_to_tskit"
missing_bin="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/missings"

declare -A names_and_muts=(
    ["Ne2e+04_r1e-08_100_1e+06"]="1e-8"
    ["Ne2e+04_r2e-08_100_1e+06"]="2e-8"
    ["Ne2e+04_r5e-08_100_1e+06"]="5e-8"
    ["Ne2e+04_r5e-09_100_1e+06"]="5e-9"
)

# Build flat list of all missing replicates
flat_list="/fs/cbsubscb18/storage/Jinmin/slurm/all_missing_list.txt"
> "$flat_list"
for name in "${!names_and_muts[@]}"; do
    missing_file="${missing_bin}/${name}.txt"
    while read -r i; do
        [ -z "$i" ] && continue
        echo "${name} ${i}" >> "$flat_list"
    done < "$missing_file"
done

# Number of array tasks and batch size
total=$(wc -l < "$flat_list")
batch_size=$(( (total + 19) / 20 ))  # 20 array tasks
start=$(( SLURM_ARRAY_TASK_ID * batch_size + 1 ))
end=$(( start + batch_size - 1 ))
[ $end -gt $total ] && end=$total

echo "Array task ${SLURM_ARRAY_TASK_ID} processing lines $start-$end of $total"
for idx in $(seq $start $end); do
    line=$(sed -n "${idx}p" "$flat_list")
    [ -z "$line" ] && continue   # skip if line is empty

    name=$(echo $line | awk '{print $1}')
    i=$(echo $line | awk '{print $2}')
    rate=${names_and_muts[$name]}

    echo "Running replicate $i for $name (mutation rate $rate) on node $(hostname)"

    vcf_in="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}.vcf.gz"
    vcf_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}.vcf"
    singer_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/singer_temp/${name}/${i}/${name}_${i}"
    tskit_out="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/tree/${name}/${i}/singer_inferred/${name}"

    mkdir -p "$(dirname "${vcf_out}")" "$(dirname "${singer_out}")" "$(dirname "${tskit_out}")"

    gunzip -c "${vcf_in}" > "${vcf_out}"

    "${singer_bin}" \
        -Ne 2e4 \
        -m "${rate}" \
        -vcf "/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/vcf/${name}/${i}/${name}_${i}" \
        -output "${singer_out}" \
        -start 0 \
        -end 1e6 \
        -thin 50 \
        -n 50 \
        -polar 0.99

    "${convert_bin}" \
        -input "${singer_out}" \
        -output "${tskit_out}" \
        -start 30 \
        -end 50 \
        -step 1

    echo "Finished replicate $i for $name"
done


echo "Array task ${SLURM_ARRAY_TASK_ID} finished."
