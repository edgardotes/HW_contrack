#!/bin/bash

#set -euo pipefail

savepath=/scratch2/edolores/era5/Z500_era5
mkdir -p ${savepath};


for year in $(seq 1990 2021); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean day
  loadpath=/mnt/climstor/ecmwf/era5/raw/PL/data/${year}/an_pl_ERA5_${year}*
  out="${savepath}/z500_day_${year}_1p0deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
    base="$(basename "${fn%.nc}")"
    output="${savepath}/${base}_z500_day_1p0deg.nc"   # distinct processed filename
    if [ ! -f "$output" ]; then
        echo $output
        cdo -P 8 -L -b F32  -remapbil,gridfile_1x1.txt  -select,name=z,level=500 $fn $output
   fi
  done

# merge only the processed files
  out="${savepath}/z500_day_${year}_1p0deg.nc"
  cdo -mergetime "${savepath}"/an_pl_ERA5_${year}*.nc "$out"
#

# optional: 
  rm -r $savepath/an_pl_ERA5_${year}*.nc

done
