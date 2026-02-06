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

savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/U500
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/2001/an_pl_ERA5_2001-08-28.nc > gridfile.txt

for year in $(seq 2020 2025); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean day
  loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_u500_monmean_${year}*
  out="${savepath}/u500_mon_${year}_0p5deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
    base="$(basename "${fn%.nc}")"
    output="${savepath}/${base}_u500_day_0p5deg.nc"   # distinct processed filename
    if [ ! -f "$output" ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt  -select,name=u $fn $output
   fi
  done

# merge only the processed files
  cdo -mergetime "${savepath}"/era5_u500_monmean_${year}*.nc "$out"
#

# optional: 
  rm -r $savepath/era5_u500_monmean_${year}*.nc

done
