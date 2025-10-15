#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=03:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

year=2024
##==extract t2m and compute mean day
#loadpath=/mnt/climstor/ecmwf/era5/raw/SFC/data/${year}/an_sfc_ERA5_${year}*
loadpath=/capstor/scratch/cscs/edolores/OBS/ERA5/tmp/era5_t2m_${year}*
#savepath=/scratch2/edolores/era5/TS_era5
savepath=/capstor/scratch/cscs/edolores/OBS/ERA5/TS_era5
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes /mnt/climstor/ecmwf/era5/raw/SFC/data/1959/an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes /mnt/climstor/ecmwf/era5/raw/PL/data/2001/an_pl_ERA5_2001-08-28.nc > gridfile.txt

out="${savepath}/t2m_day_${year}_0p5deg.nc"

###
for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymean -select,name=t2m $fn $output
   fi
done

cdo -mergetime $savepath/era5_t2m_${year}*.nc $out

# optional: 
rm -r $savepath/era5_t2m_${year}*.nc

