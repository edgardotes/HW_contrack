#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=08:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

###########
# ERA5
#############

savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/Z500_era5
mkdir -p ${savepath};

for year in $(seq 2009 2019); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean day
  loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_z500_${year}*
  out="${savepath}/z500_day_${year}_0p5deg.nc" 

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
    base="$(basename "${fn%.nc}")"
    output="${savepath}/${base}_z500_day_0p5deg.nc"   # distinct processed filename
    if [ ! -f "$output" ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt   -daymean -select,name=z,level=500 $fn $output
   fi
  done

# merge only the processed files
  cdo -mergetime "${savepath}"/era5_z500_${year}*.nc $out

# optional: 
  rm -r $savepath/era5_z500_${year}*.nc

done
