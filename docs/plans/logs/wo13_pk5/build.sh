#!/bin/bash
# wo13-pk5: Step-1's M6 rebuilt in a scratch root from ~/ctx_tier1's gather and DE outputs.
#   bash build.sh <scratch root>
# The root links the real data, resources, 01_gather and 05_de (read only by ctx function)
# and holds its own 06_function.
real=/Users/rytis/ctx_tier1
here=$(cd "$(dirname "$0")" && pwd)
root=$1
mkdir -p "$root/outputs"
for d in data resources; do [ -e "$root/$d" ] || ln -s "$real/$d" "$root/$d"; done
for d in 01_gather 02_qc 05_de; do
  [ -e "$root/outputs/$d" ] || ln -s "$real/outputs/$d" "$root/outputs/$d"
done
/usr/bin/time -l "$here/py.sh" python -m ctx --root "$root" function --dataset tier1 \
  --dataset-file conf/datasets/tier1.toml
