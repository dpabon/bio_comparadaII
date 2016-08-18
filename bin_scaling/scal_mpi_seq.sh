#!/bin/bash
cd /home/dane/MEGAsync/bio_comparadaII/scaling/
mpirun -np 4 rb-mpi scal_mpi_mcmc.rev
rb scal_mpi_map.rev
