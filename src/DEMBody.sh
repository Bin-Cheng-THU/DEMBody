#!/bin/bash
ulimit -s unlimited
export OMP_STACKSIZE=1000m
srun -N 1 -n 1 -c 24 DEMBody