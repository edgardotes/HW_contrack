#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

###########
# ERA5
#############

year=2024

##==extract t2m and compute mean day
#loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_t850_${year}*
loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_z500_${year}*

#savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/T850_era5
savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/Z500_era5
mkdir -p ${savepath};

#out="${savepath}/t850_day_${year}_0p5deg.nc"
out="${savepath}/z500_day_${year}_0p5deg.nc"

###
for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
#        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymean -select,name=t,level=850 $fn $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymean -select,name=z,level=500 $fn $output
   fi
done

#cdo -mergetime $savepath/era5_t850_${year}*.nc $out
cdo -mergetime $savepath/era5_z500_${year}*.nc $out


# optional: 
#rm -r $savepath/era5_t850_${year}*.nc
rm -r $savepath/era5_z500_${year}*.nc

