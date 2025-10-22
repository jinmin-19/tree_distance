#!/bin/bash

# List of scenario names
names=("Ne2e+04_r1e-08_100_1e+06" "Ne2e+04_r2e-08_100_1e+06" "Ne2e+04_r5e-08_100_1e+06" "Ne2e+04_r5e-09_100_1e+06")

src_base="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/diff_rate/tree"
dst_base="/fs/cbsubscb18/storage/Jinmin/tree_dis/simulations/small_diff_rate"

for name in "${names[@]}"; do
    for i in {0..99}; do
        src="$src_base/$name/$i"
        dst="$dst_base/$name/$i"

        mkdir -p "$dst"

        # Copy all files in the folder except singer_inferred
        find "$src" -maxdepth 1 -type f -exec cp {} "$dst/" \;

        # Copy only *49.trees from singer_inferred
        mkdir -p "$dst/singer_inferred"
        cp "$src/singer_inferred/"*49.trees "$dst/singer_inferred/" 2>/dev/null
    done
done
