#!/bin/bash
rm pastis-*.log ; rm slurm.ipu* ; cd build/ && make -j$(nproc) && cd .. && sbatch runpastis
