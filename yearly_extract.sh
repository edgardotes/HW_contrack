#!/bin/bash

#set -euo pipefail

#year=1963
savepath=/scratch2/edolores/era5/TS_era5
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/2001/an_pl_ERA5_2001-08-28.nc > gridfile.txt

for year in $(seq 1964 1970); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean day
  loadpath=/mnt/climstor/ecmwf/era5/raw/SFC/data/${year}/an_sfc_ERA5*
  out="${savepath}/t2m_day_${year}_0p5deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymean -select,name=t2m $fn $output
   fi
  done

  cdo -mergetime $savepath/an_sfc_ERA5_${year}*.nc $out


# optional: 
  rm -r $savepath/an_sfc_ERA5_${year}*.nc

done
