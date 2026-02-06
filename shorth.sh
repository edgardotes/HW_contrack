#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=04:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

# Keep only essentials; adjust variable names to your file (common ICON names shown)
#ncks -O \
#  -v clon,clat,clon_vertices,clat_vertices,cell_area \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/icon_grid_0056_R02B10_G.nc \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/grid_R02B10_light.nc

# Turn on netCDF4 compression and reasonable chunking
#nccopy -k netCDF-4 -d 3 -s \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/grid_R02B10_light.nc \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/grid_R02B10_light_c.nc

### STEP 2
### For 50 km (0.5°):
#cdo -P 8 gencon,r720x360 \
#  -setgrid,/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/grid_R02B10_light_c.nc \
#  /capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_2/tot_prec_20210602T000000Z.nc \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/weights_R02B10_to_0p5deg.nc

# For 12.5 km (0.125°):
#cdo -P 8 gencon,r2880x1440 \
#  -setgrid,/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/grid_R02B10_light_c.nc \
#  /capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_2/tot_prec_20210602T000000Z.nc \
#  /capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/weights_R02B10_to_0p125deg.nc

### STEP 3

export OMP_NUM_THREADS=8   # keep threads moderate to avoid RAM spikes

var=tot_prec
FILE_PATH=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_2/$var
SAVE_PATH=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM
#WEIGHTS=/capstor/scratch/cscs/edolores/old_icon/weights_R02B10_to_0p125deg.nc  
WEIGHTS=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/weights_R02B10_to_0p5deg.nc


echo "[$WEIGHTS]"
ls -lh "$WEIGHTS"
cdo -showname -remap,"$WEIGHTS" -selname,tot_prec -seltimestep,1 \
    -setgrid,${SAVE_PATH}/grid_R02B10_light_c.nc \
    /capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_2/tot_prec_20200120T000000Z.nc \
    /dev/null

