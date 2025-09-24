#!/bin/bash

#set -euo pipefail

year=1959
##==extract t2m and compute mean day
loadpath=/mnt/climstor/ecmwf/era5/raw/PL/data/${year}/an_pl_ERA5_${year}*
savepath=/scratch2/edolores/era5/T850_era5
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/2001/an_pl_ERA5_2001-08-28.nc > gridfile.txt

out="${savepath}/t850_day_${year}_0p5deg.nc"


###
for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -daymean -select,name=t,level=850 $fn $output
   fi
done

cdo -mergetime $savepath/an_pl_ERA5_${year}*.nc $out


# optional: 
rm -r $savepath/an_pl_ERA5_${year}*.nc

