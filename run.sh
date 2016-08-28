#!/bin/bash

DIR=../

./build.sh

for i in $(seq 9199); do
    idx=$(printf %.5d $i)
    ./find_oscillons $DIR/cuda/output-20160817091803/rho_${idx}.bin /tmp/out$1 128 8 8 1 /tmp/stats$1 1>/dev/null
    ./extract_profile $DIR/cuda/output-20160817091803/phi_${idx}.bin $1 /tmp/out$1 > /tmp/prof$1
    nrlines=$(wc -l /tmp/prof$1 | cut -f1 -d' ')
    head -n $((nrlines / 2)) /tmp/prof$1 | tail -n 1
done
