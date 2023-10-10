#!/bin/bash
#SBATCH --job-name=build_install_model
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trobinet@stanford.edu
#SBATCH --nodes=5
#SBATCH --tasks=1
#SBATCH --time=03:00:00
#SBATCH --output=/shared/models/GEOSldas_pso_g1_ksat/GEOSldas/build_install_logs/sbatch_output_%j.log
#SBATCH --partition=catch-m6a4xl-spot

echo "running!"

# run the code
./build_install_model.sh

echo "done!"
