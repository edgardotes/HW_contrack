from pyproj import Proj, Transformer, transform, CRS
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
#import xesmf as xe
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import glob
import os
import pandas as pd
import time
from itertools import groupby
from scipy.stats import linregress
from scipy import ndimage
from contrack import contrack
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely
import shapely.ops as sops
import shapely.vectorized as sv

# ----------------------------
# Helpers
# ----------------------------
# --- keep this if you already have it ---
def drop_feb29(da: xr.DataArray) -> xr.DataArray:
    is_feb29 = (da['time'].dt.month == 2) & (da['time'].dt.day == 29)
    return da.sel(time=~is_feb29)

# NEW: day-of-year remapped to 1..365 (no leap 29 Feb)
def dayofyear_365(time_index):
    doy = time_index.dt.dayofyear
    leap = time_index.dt.is_leap_year
    # shift all days after Feb 28 in leap years by -1
    return xr.where(leap & (doy > 59), doy - 1, doy)

def circular_day_distance(doy, center, period=365):
    return np.abs(((doy - center + period//2) % period) - period//2)

def compute_tr90d(
    tmax: xr.DataArray,
    ref_start: str = "1961-01-01",
    ref_end: str = "1990-12-31",
    window_halfwidth: int = 15,
    q: float = 0.90,  # fraction
) -> xr.DataArray:
    # remove Feb 29 and slice reference
    tmax_noleap = drop_feb29(tmax)
    ref = drop_feb29(tmax.sel(time=slice(ref_start, ref_end)))

    # base chunking (tune to your system)
    if ref.chunks is None:
        ref = ref.chunk({"time": 365, "lat": 90, "lon": 180})

    # IMPORTANT: use no-leap day-of-year for window selection (1..365)
    ref_doy365 = dayofyear_365(ref["time"])

    doy_vals = np.arange(1, 366)
    out = []
    for k in doy_vals:
        dist = circular_day_distance(ref_doy365, k, period=365)
        sel = ref.where(dist <= window_halfwidth, drop=True)
        # single-chunk small window along time so quantile is happy
        sel = sel.chunk({"time": -1})
        tr_k = sel.quantile(q, dim="time", skipna=True)
        out.append(tr_k)

    tr = xr.concat(out, dim=xr.DataArray(doy_vals, dims="dayofyear", name="dayofyear"))
    tr.name = f"Tr{int(q*100)}d"
    tr.attrs.update({
        "description": f"{int(q*100)}th percentile of daily Tmax (±{window_halfwidth} days) over {ref_start[:4]}–{ref_end[:4]}",
        "units": tmax.attrs.get("units", ""),
    })
    return tr

def label_runs_1d(a):
    labels, _ = ndimage.label(a.astype(np.uint8))
    return labels

def run_lengths_from_labels_1d(labels):
    if labels.size == 0:
        return labels
    max_lab = labels.max()
    if max_lab == 0:
        return np.zeros_like(labels)
    counts = np.bincount(labels, minlength=max_lab + 1)
    lengths = counts[labels]
    lengths[labels == 0] = 0
    return lengths

def detect_heatwaves(tmax: xr.DataArray, tr_by_doy: xr.DataArray, min_duration: int = 3):
    assert min_duration == 3, "This fast path is coded for min_duration=3. Ask me if you need a general N."

    # work with no-leap series
    t_noleap = drop_feb29(tmax)

    # 1..365 mapping (no 366s)
    doy365 = dayofyear_365(t_noleap.time)

    # expand Tr to full time (dask-friendly)
    tr_full = tr_by_doy.sel(dayofyear=doy365).transpose("time", "lat", "lon")
    tr_full.name = tr_by_doy.name

    # --- NEW: ensure both operands have large-enough time chunks [error of 14.10.2025]
    # (vectorized selection above often produces time chunks of size 1)
    T = t_noleap.sizes["time"]
    t_noleap = t_noleap.chunk({"time": -1})                 # or e.g. 365
    tr_full  = tr_full.chunk({"time": -1})                  # match or merge-friendly

    # hot days (exceed threshold)
    hot_mask = (t_noleap > tr_full).rename("hot_day")

    # ---- NEW: rolling trick, avoids apply_ufunc & single-chunk time ----
    end3 = hot_mask.rolling(time=3, min_periods=3).sum() == 3
    hw_mask = end3 | end3.shift(time=1, fill_value=False) | end3.shift(time=2, fill_value=False)
    hw_mask = hw_mask.rename("heatwave_day")

    return hw_mask, hot_mask, tr_full

# ----------------------------
# Seasonal split of heatwave events
# ----------------------------
def _season_masks_from_profile(hw_1d, month_1d):
    """
    Per-gridpoint numpy routine (time only).
    Inputs:
      hw_1d:   (T,) bool — heatwave day mask (already min_duration>=3)
      month_1d:(T,) int  — calendar month 1..12
    Returns:
      out: (4, T) bool for seasons [DJF, MAM, JJA, SON]
    """
    # seasons definition
    core = [
        {12, 1, 2},      # DJF
        {3, 4, 5},       # MAM
        {6, 7, 8},       # JJA
        {9, 10, 11},     # SON
    ]
    window = [
        {11, 12, 1, 2, 3},  # DJF window: Nov–Mar
        {2, 3, 4, 5, 6},    # MAM window: Feb–Jun
        {5, 6, 7, 8, 9},    # JJA window: May–Sep
        {8, 9, 10, 11, 12}, # SON window: Aug–Dec
    ]

    out = np.zeros((4, hw_1d.size), dtype=bool)
    if hw_1d.size == 0 or not np.any(hw_1d):
        return out

    # label contiguous heatwave runs
    labels, nlab = ndimage.label(hw_1d.astype(np.uint8))
    if nlab == 0:
        return out

    # iterate events
    for lab in range(1, nlab + 1):
        idx = np.nonzero(labels == lab)[0]
        if idx.size == 0:
            continue

        # count overlap in core months for each season
        counts_core = []
        for s in range(4):
            counts_core.append(np.isin(month_1d[idx], list(core[s])).sum())
        counts_core = np.asarray(counts_core)

        # choose season: argmax over core counts (ties -> first in [DJF,MAM,JJA,SON])
        s_star = int(np.argmax(counts_core))

        # keep only the days of this event that lie in the 5-month window of the assigned season
        keep = np.isin(month_1d[idx], list(window[s_star]))
        if np.any(keep):
            out[s_star, idx[keep]] = True

    return out

def split_heatwaves_by_season(hw_mask: xr.DataArray) -> xr.DataArray:
    hw = drop_feb29(hw_mask)
    months = hw["time"].dt.month

    # Keep your spatial chunking; allow the gufunc to rechunk time as needed
    if hw.chunks is None:
        hw = hw.chunk({"time": 365, "lat": 90, "lon": 180})

    seasons = np.array(["DJF", "MAM", "JJA", "SON"])

    season_hw = xr.apply_ufunc(
        _season_masks_from_profile,
        hw.transpose("time", "lat", "lon"),
        months,
        input_core_dims=[["time"], ["time"]],
        output_core_dims=[["season", "time"]],
        output_dtypes=[bool],
        dask="parallelized",
        vectorize=True,
        output_sizes={"season": 4},
        dask_gufunc_kwargs={"allow_rechunk": True},  # <- key change
    )

    season_hw = season_hw.assign_coords(season=("season", seasons))
    season_hw = season_hw.transpose("season", "time", "lat", "lon")
    season_hw.name = "heatwave_day_by_season"
    season_hw.attrs.update({
        "description": (
            "Heatwave days split by season using 5-month windows around each 3-month season. "
            "Events are assigned to the season where they have the most days in the core months; "
            "for that season we keep only the event's days within that season's 5-month window."
        )
    })
    return season_hw

###opening saved data
hw_mask=xr.open_dataset('/capstor/scratch/cscs/edolores/OBS/ERA5/Heatwaves/Russo/hw_mask.nc')['heatwave_day']
hot_mask=xr.open_dataset('/capstor/scratch/cscs/edolores/OBS/ERA5/Heatwaves/Russo/hot_mask.nc')['hot_day']
tr_full=xr.open_dataset('/capstor/scratch/cscs/edolores/OBS/ERA5/Heatwaves/Russo/tr_full.nc')['Tr90d']

# 3) Season-split - returns a [season, time, lat, lon] boolean mask
season_hw = split_heatwaves_by_season(hw_mask)

# Example: counts full period
### average per season
jja = season_hw.sel(season="JJA")                 # includes May-Sep by design
jja_yearly = jja.groupby("time.year").sum("time") # days per year at each gridpoint
jja_yearly_mean = jja_yearly.mean("year").compute()         # average per year
jja_yearly_mean.attrs["units"] = "days per season"
jja_yearly_mean.attrs["long_name"] = "Mean JJA heatwave days per season"
print(jja_yearly_mean.max().values)

jja_yearly_mean.to_netcdf('/capstor/scratch/cscs/edolores/OBS/ERA5/Heatwaves/Russo/jja_yearly_mean.nc')
