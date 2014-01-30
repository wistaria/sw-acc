#!/bin/sh

rm -f omp.dat
for nproc in `seq 1 1 15`; do
  OMP_NUM_THREADS=$nproc ./rand_pi_omp 100000000 10 >> omp.dat
done
