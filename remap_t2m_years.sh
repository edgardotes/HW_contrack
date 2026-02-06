#!/bin/bash
#
#SBATCH --account=cwp03
#SBATCH --job-name=test
#SBATCH --time=02:00:00
#SBATCH --nodes=4
#SBATCH --output=test.o
#SBATCH --error=test.e

## tmax, tmin, t2m (tmean)
var=t2m   # choose: tmin | tmax | t2m

###ICON 2.5 km
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed/out_7/
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/TS

###ICON 10 km
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08/out_7/
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM/TS
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM/TMAX
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM/TMIN

###ICON 10 km ON
loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08_shallow_only/out_7/
savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM-ON/TS
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM-ON/TMAX
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_10KM-ON/TMIN

###ICON 40 km
#loadpath=/capstor/store1/cscs/userlab/cwp03/zemanc/Data_Dyamond_PostProcessed_R02B08_shallow_only/out_7/
#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_40KM/TS

# choose CDO operator based on var
case "$var" in
  tmin)
    cdo_op=daymin
    ;;
  tmax)
    cdo_op=daymax
    ;;
  t2m)
    cdo_op=daymean
    ;;
  *)
    echo "Unknown var='$var'. Use one of: tmin, tmax, t2m."
    exit 1
    ;;
esac

#savepath=/capstor/scratch/cscs/edolores/EXCLAIM/ICON_2KM/$(echo "$var" | tr '[:lower:]' '[:upper:]')
mkdir -p "${savepath}"

for year in $(seq 2020 2024); do
  echo "=== Processing $year ($var / $cdo_op) ==="

  # input files (still t_2m in your filenames)
  loadpath_year=${loadpath}/remap_t_2m_${year}*
  out="${savepath}/${var}_day_${year}_0p5deg.nc"

  # skip if final yearly file exists
  if [ -f "$out" ]; then
    echo "Already have $out, skipping"
    continue
  fi

  # process each file: subset + remap + daily stat
  for fn in $loadpath_year; do
    base="$(basename "${fn%.nc}")"
    # give output files a var-specific name
    output="${savepath}/${base}_${var}_day_0p5deg.nc"

    if [ ! -f "$output" ]; then
      echo "  -> creating $output"
      cdo -P 8 -L -b F32 \
        -sellonlatbox,-180,180,-90,90 \
        -remapbil,gridfile.txt \
        -${cdo_op} \
        -select,name=t_2m \
        "$fn" "$output"
    fi
  done

  # merge only the processed files for that year & var
  cdo -mergetime "${savepath}/remap_t_2m_${year}"*"_${var}_day_0p5deg.nc" "$out"

  # optional: clean up intermediates
  rm -f "${savepath}/remap_t_2m_${year}"*"_${var}_day_0p5deg.nc"

done
