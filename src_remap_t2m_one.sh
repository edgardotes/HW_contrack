#!/bin/bash

year=2022
##==extract t2m and compute mean day
###ICON 2.5 km
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_7/remap_t_2m_${year}*
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/TS
###ICON 10 km
loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08/out_7/remap_t_2m_${year}*
savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM/TS
###ICON 10 km ON
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08_shallow_only/out_7/remap_t_2m_${year}*
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM-ON/TS
###ICON 40 km
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08_shallow_only/out_7/remap_t_2m_${year}*
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_40KM/TS

mkdir -p ${savepath};

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

