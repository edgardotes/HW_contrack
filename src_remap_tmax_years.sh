#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=02:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

###create grid
##0.25
#cdo griddes an_sfc_ERA5_1959-08-26.nc > gridfile.txt
###0.50
#cdo griddes an_pl_ERA5_2001-08-28.nc > gridfile.txt

## tmax, tmin, t2m (tmean)
var=tmax
savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/TMAX
mkdir -p ${savepath};

for year in $(seq 2020 2024); do
  echo "=== Processing $year ==="

##==extract t2m and compute mean/max/min day
  loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_7/remap_t_2m_${year}*
#  out="${savepath}/t2m_day_${year}_0p5deg.nc"
  out="${savepath}/${var}_day_${year}_0p5deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

###
  for fn in $loadpath; do
    base="$(basename "${fn%.nc}")"
    output="${savepath}/${base}_t2m_day_0p5deg.nc"   # distinct processed filename
    if [ ! -f "$output" ]; then
        echo $output
        cdo -P 8 -L -b F32 -sellonlatbox,-180,180,0,90  -remapbil,gridfile.txt -daymax -select,name=t_2m $fn $output

   fi
  done

# merge only the processed files
  cdo -mergetime "${savepath}"/remap_t_2m_${year}*.nc "$out"
#

# optional: 
  rm -r $savepath/remap_t_2m_${year}*.nc

done
