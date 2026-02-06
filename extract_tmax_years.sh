#!/bin/bash

###########
# ERA5
#############

#year=1963
savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/TMAX_era5
mkdir -p ${savepath};

for year in $(seq 2024 2025); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean day
#  loadpath=/mnt/climstor/ecmwf/era5/raw/SFC/data/${year}/an_sfc_ERA5_${year}*
  loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_t2m_${year}*
#  out="${savepath}/t2m_day_${year}_0p5deg.nc"
  out="${savepath}/tmax_day_${year}_0p5deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
#    fn_new=$savepath/$(basename $fn)
#    output=${fn_new}
    base="$(basename "${fn%.nc}")"
    output="${savepath}/${base}_t2m_day_0p5deg.nc"   # distinct processed filename
    if [ ! -f "$output" ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymax -select,name=t2m $fn $output
   fi
  done

# merge only the processed files
  cdo -mergetime "${savepath}"/era5_t2m_${year}*.nc "$out"
#

# optional: 
  rm -r $savepath/era5_t2m_${year}*.nc

done
