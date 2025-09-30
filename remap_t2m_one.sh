#!/bin/bash

year=2023
##==extract t2m and compute mean day
loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_7/remap_t_2m_${year}*
savepath=/capstor/scratch/cscs/edolores/ICON_2KM/TS_era5
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes an_pl_ERA5_2001-08-28.nc > gridfile.txt

out="${savepath}/t2m_day_${year}_0p5deg.nc"


###
for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymean -select,name=t_2m $fn $output
   fi
done

cdo -mergetime $savepath/remap_t_2m_${year}*.nc $out


# optional: 
rm -r $savepath/remap_t_2m_${year}*.nc

