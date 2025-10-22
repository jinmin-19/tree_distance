#!/bin/bash

names=(
  "Ne2e+04_r1e-08_100_1e+06"
  "Ne2e+04_r2e-08_100_1e+06"
  "Ne2e+04_r5e-08_100_1e+06"
  "Ne2e+04_r5e-09_100_1e+06"
)

for name in "${names[@]}"; do
  base="/Users/jinminli/research/tree_distance/simulations/diff_rates/${name}"
  outlist="/Users/jinminli/research/tree_distance/simulations/diff_rates/missings/${name}.txt"

  mkdir -p "$(dirname "$outlist")"
  > "$outlist"

  for i in {0..99}; do
    if [ ! -f "${base}/${i}/singer_inferred/${name}_49.trees" ]; then
      echo "$i" >> "$outlist"
    fi
  done

  echo "Missing replicates listed in $outlist"
done
