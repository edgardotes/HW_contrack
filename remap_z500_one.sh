#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

year=2024
##==extract Z from remaped, 3 hourly with 37 levels
loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_1_1/remap_geopot_${year}*
###UWIND
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_1_4/remap_u_${year}*
###Temperature
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_1_3/remap_temp_${year}*

### saving
savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/Z500
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/U500
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/T850
mkdir -p ${savepath};

###create grid
##0.25
#cdo griddes an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes an_pl_ERA5_2001-08-28.nc > gridfile.txt

### file name output
out="${savepath}/z500_day_${year}_0p5deg.nc"
#out="${savepath}/u500_day_${year}_0p5deg.nc"
#out="${savepath}/t850_day_${year}_0p5deg.nc"

###
for fn in $loadpath; do
    fn_new=$savepath/$(basename $fn)
    output=${fn_new}
    if [ ! -f ${output} ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,-90,90  -remapbil,gridfile.txt -daymean -select,name=geopot,level=50000 $fn $output
#        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,-90,90  -remapbil,gridfile.txt -daymean -select,name=u,level=50000 $fn $output
#        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,-90,90  -remapbil,gridfile.txt -daymean -select,name=temp,level=85000 $fn $output
   fi
done

cdo -mergetime $savepath/remap_geopot_${year}*.nc $out
#cdo -mergetime $savepath/remap_u_${year}*.nc $out
#cdo -mergetime $savepath/remap_temp_${year}*.nc $out


# optional: 
rm -r $savepath/remap_geopot_${year}*.nc
#rm -r $savepath/remap_u_${year}*.nc
#rm -r $savepath/remap_temp_${year}*.nc

