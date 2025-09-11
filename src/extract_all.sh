#!/bin/bash


##==extract z500 and compute mean day
loadpath=/mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5*
savepath=/home/edgar/NEXTGEMS/data/TS_era5

###create grid
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/an_pl_ERA5_2001-08-28.nc > gridfile.txt

#### create folder to save data
if [ ! -d ${savepath}  ]; then
   mkdir -p ${savepath};
fi

for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f output ]; then
        echo $output
#        cdo  -daymean -select,name=t,level=500 $fn $output
        cdo  remapbil,gridfile.txt -daymean -select,name=t2m $fn $output
   fi
done


