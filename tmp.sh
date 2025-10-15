#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=04:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

#uenv start netcdf-tools/2024:v1 --view=netcdf-tools:default
source ~/.env_icon/bin/activate

### download era5 data
#python get_era5_t2m.py

### compute seasonal freqeuncy of heatwaves
python get_hw_season.py
