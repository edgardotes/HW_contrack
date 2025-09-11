#!/bin/bash

###
#yy0=1970
#yyn=2019

##==extract z500 and compute mean day
#loadpath=/mnt/climstor/ecmwf/era5/raw/PL/data/an_pl_ERA5*
loadpath=/mnt/climstor/ecmwf/era5/raw/PL/data/an_pl_ERA5_197*
savepath=/home/edolores/NEXTGEMS/tools/weather_regimes/7YWR/data/z500
outfile=era5_z500_natl_daily_1970-1979_2deg.nc

## sellonlatbox: sellonlatbox,lon1,lon2,lat1,lat2
## NA region                  -85.5, 45, 27, 81 

###create grid
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
#cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/an_pl_ERA5_2001-08-28.nc > gridfile.txt
#gfortran -c griddes gridfile.f  ### change some parameters!!!
#./griddes

#### create folder to save data
if [ ! -d ${savepath}  ]; then
   mkdir -p ${savepath};
fi

##cdo mergetime -sellonlatbox,-85.5,45,27,81 -remapbil,gridfile.txt -daymean -select,name=z,level=500 $loadpath $outfile
cdo mergetime -sellonlatbox,-80,40,30,90 -remapbil,gridfile.txt -daymean -select,name=z,level=500 $loadpath $outfile
