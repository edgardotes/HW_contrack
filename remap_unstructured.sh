#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=04:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

var=t_2m #tmax_2m
GRID_PATH=/capstor/scratch/cscs/edolores/old_icon
FILE_PATH=/capstor/scratch/cscs/ppothapa/for_edgar/$var
SAVE_PATH=/capstor/scratch/cscs/edolores/old_icon/$var

mkdir -p ${SAVE_PATH};

## max value
#for year in $(seq 2006 2017); do
#  echo "=== Processing $year ==="
#  cdo -P 8 -f nc remapcon,r720x360 \
#  -daymax \
#  -setgrid,$GRID_PATH/icon_grid_0052_R02B06_G.nc \
#  $FILE_PATH/${var}_${year}_compressed.nc \
#  $SAVE_PATH/remap_${var}_${year}_0.5deg.nc
#
#done
## mean value
for year in $(seq 2006 2017); do
  echo "=== Processing $year ==="
  cdo -P 8 -f nc remapcon,r720x360 \
  -daymean \
  -setgrid,$GRID_PATH/icon_grid_0052_R02B06_G.nc \
  $FILE_PATH/${var}_${year}_compressed.nc \
  $SAVE_PATH/remap_${var}_${year}_0.5deg.nc

done
