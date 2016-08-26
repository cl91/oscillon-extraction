#!/bin/bash

./build.sh

for i in $(seq 9199); do
    idx=$(printf %.5d $i)
    ./find_oscillons /mnt/nesi/cuda/output-20160817091803/rho_${idx}.bin /tmp/out 128 8 8 1 /tmp/stats 1>/dev/null
    ./extract_profile /mnt/nesi/cuda/output-20160817091803/phi_${idx}.bin 2 /tmp/out > /tmp/prof
    nrlines=$(wc -l /tmp/prof | cut -f1 -d' ')
    head -n $((nrlines / 2)) /tmp/prof | tail -n 1
done
