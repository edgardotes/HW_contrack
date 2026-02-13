# ----------------------------
# Helpers
# ----------------------------
# Traditional grid cell HW identification
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
import matplotlib.colors as colors
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

        # --- NEW: ensure both operands have large-enough time chunks [error of 10.10.2025]
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
    """
    Assign each heatwave event to a season as per the 5-month window rule
    and return a per-season boolean mask of heatwave days.

    Parameters
    ----------
    hw_mask : xr.DataArray [time, lat, lon], bool
        Output `heatwave_day` from detect_heatwaves (already no-leap).

    Returns
    -------
    season_hw : xr.DataArray [season, time, lat, lon], bool
        Per-season masks with seasons = ['DJF','MAM','JJA','SON'].
        No event is double-counted across seasons.
    """
    # ensure Feb 29 removed (your detect_heatwaves already works no-leap, but be safe)
    hw = drop_feb29(hw_mask)

    # month index (1..12) aligned with hw time
    months = hw["time"].dt.month

    # make sure chunking is friendly (time-wise vector ops)
    if hw.chunks is None:
        hw = hw.chunk({"time": 365, "lat": 90, "lon": 180})

    seasons = np.array(["DJF", "MAM", "JJA", "SON"])

    # apply per pixel along time
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
    )

    # tidy up coords and attrs
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

# ===============================
# Option 2: Full-event durations per season-year
# ===============================
def _events_full_len_assign_season(hw_1d, month_1d, year_1d):
    """
    Numpy routine (time-only) to:
      - label contiguous heatwave runs,
      - compute each event's FULL duration (end - start + 1),
      - assign the WHOLE event to exactly one season by core-month majority,
      - define a season-year for the event,
      - emit per-time arrays that are ONLY filled at event END days:
          * length_at_end[T] : int (duration days), else 0
          * season_idx_at_end[T] : {0..3} for [DJF,MAM,JJA,SON], else -1
          * syear_at_end[T] : int season-year, else -1
    """
    T = hw_1d.size
    if T == 0 or not np.any(hw_1d):
        return (np.zeros(T, dtype=np.int16),
                -np.ones(T, dtype=np.int8),
                -np.ones(T, dtype=np.int32))

    # seasons & core-months
    seasons = ["DJF", "MAM", "JJA", "SON"]
    core = [
        {12, 1, 2},      # DJF
        {3, 4, 5},       # MAM
        {6, 7, 8},       # JJA
        {9, 10, 11},     # SON
    ]

    # label runs
    labels, nlab = ndimage.label(hw_1d.astype(np.uint8))

    length_at_end = np.zeros(T, dtype=np.int16)
    season_idx_at_end = -np.ones(T, dtype=np.int8)
    syear_at_end = -np.ones(T, dtype=np.int32)

    for lab in range(1, nlab + 1):
        idx = np.nonzero(labels == lab)[0]
        if idx.size == 0:
            continue

        # FULL duration (days)
        dur = idx.size

        # assign season by core-month majority (ties -> DJF,MAM,JJA,SON order)
        months_evt = month_1d[idx]
        counts = np.array([np.isin(months_evt, list(core[s])).sum() for s in range(4)])
        s_star = int(np.argmax(counts))  # 0..3

        # season-year via event midpoint (robust for cross-year events)
        mid = idx[len(idx) // 2]
        y_mid = int(year_1d[mid])
        m_mid = int(month_1d[mid])
        if s_star == 0 and m_mid == 12:  # DJF & midpoint in Dec -> next year
            syear = y_mid + 1
        else:
            syear = y_mid

        # write only at the event END day
        end_day = idx[-1]
        length_at_end[end_day] = dur
        season_idx_at_end[end_day] = s_star
        syear_at_end[end_day] = syear

    return length_at_end, season_idx_at_end, syear_at_end


def durations_full_event_per_season_year(hw_mask: xr.DataArray) -> xr.Dataset:
    """
    Compute full-event heatwave durations and aggregate by (season, season-year).

    Parameters
    ----------
    hw_mask : xr.DataArray [time, lat, lon], bool
        Heatwave day mask with min_duration already imposed (e.g., from detect_heatwaves).
        Time should be no-leap (Feb 29 removed). If not, we remove it here.

    Returns
    -------
    ds : xr.Dataset with variables:
        - mean_duration[season, season_year, lat, lon] : days per event
        - max_duration [season, season_year, lat, lon] : days
        - event_count  [season, season_year, lat, lon] : events per season-year
        - total_days   [season, season_year, lat, lon] : sum of full-event durations
      seasons = ['DJF','MAM','JJA','SON'].
      Units set in attrs.
    """
    # ensure no-leap and friendly chunking
    x = drop_feb29(hw_mask).astype(bool)
    if x.chunks is None:
        x = x.chunk({"time": 365, "lat": 90, "lon": 180})

    months = x["time"].dt.month.astype("int16")
    years  = x["time"].dt.year.astype("int32")

    # vectorized per-pixel apply over time
    (length_at_end, season_idx_at_end, syear_at_end) = xr.apply_ufunc(
        _events_full_len_assign_season,
        x.transpose("time", "lat", "lon"),
        months,
        years,
        input_core_dims=[["time"], ["time"], ["time"]],
        output_core_dims=[["time"], ["time"], ["time"]],
        output_dtypes=[np.int16, np.int8, np.int32],
        vectorize=True,
        dask="parallelized",
    )

    # set names & coords
    length_at_end = length_at_end.rename("length_at_end").transpose("time", "lat", "lon")
    season_idx_at_end = season_idx_at_end.rename("season_idx_at_end").transpose("time", "lat", "lon")
    syear_at_end = syear_at_end.rename("season_year_at_end").transpose("time", "lat", "lon")

    # Convenience: build a categorical-like selector per season
    seasons = np.array(["DJF","MAM","JJA","SON"])

    out = []

    for s_idx, s_name in enumerate(seasons):
        sel = (season_idx_at_end == s_idx)

        # keep only ends for this season
        len_sel = xr.where(sel, length_at_end, np.nan)

        # labels (season-year) — compute ONLY the labels to avoid the groupby-on-chunked error
        syear_labels = xr.where(sel, syear_at_end, np.nan).astype("float64").compute()

        # attach labels as a coordinate along time
#        len_lab = len_sel.assign_coords(season_year=("time", syear_labels))
#        len_lab = len_sel.assign_coords(season_year=("time", syear_labels.data))
        len_lab = len_sel.assign_coords(season_year=(len_sel.dims, syear_labels.data))

        # aggregate per season_year (xarray will create a "season_year" dimension)
        mean_dur = len_lab.groupby("season_year").mean().rename("mean_duration")

        max_dur  = len_lab.groupby("season_year").max().rename("max_duration")

        

        # event count: count of ends (True) per season_year
#        ends_lab = sel.assign_coords(season_year=("time", syear_labels))
#        ends_lab = sel.assign_coords(season_year=("time", syear_labels.data))
        ends_lab = sel.assign_coords(season_year=(sel.dims, syear_labels.data))

#        n_evt = ends_lab.groupby("season_year").sum("time").rename("event_count")
        n_evt    = ends_lab.groupby("season_year").sum().rename("event_count")
        
        # total days: sum of full-event durations at ends
#        total_days = len_lab.groupby("season_year").sum("time").rename("total_days")
        total_days = len_lab.groupby("season_year").sum().rename("total_days")

        # stack into one dataset and tag the season
        ds_s = xr.Dataset(
            dict(mean_duration=mean_dur, max_duration=max_dur,
                 event_count=n_evt, total_days=total_days)
        ).assign_coords(season=s_name).expand_dims("season")

        out.append(ds_s)

    ds = xr.concat(out, dim="season").sortby("season_year")
    ds = ds.assign_coords(season=("season", seasons))

    # tidy attrs / units
    ds["mean_duration"].attrs.update({
        "long_name": "Mean full-event duration",
        "units": "days per event",
        "description": "Average duration of heatwave events (full length), grouped by season and season-year."
    })
    ds["max_duration"].attrs.update({
        "long_name": "Max full-event duration",
        "units": "days",
    })
    ds["event_count"].attrs.update({
        "long_name": "Number of heatwave events",
        "units": "events per season-year",
    })
    ds["total_days"].attrs.update({
        "long_name": "Total heatwave days from full-event durations",
        "units": "days",
        "note": "Sum of full-event durations for events assigned to the season-year.",
    })

    return ds

def make_land_mask_like(da_2d: xr.DataArray) -> xr.DataArray:
    # discover lat/lon names
    lat_name = next(n for n in da_2d.coords if n.lower().startswith('lat'))
    lon_name = next(n for n in da_2d.coords if n.lower().startswith('lon'))

    # broadcast to 2D lon/lat with consistent dim order
    lon2d, lat2d = xr.broadcast(da_2d[lon_name], da_2d[lat_name])

    # wrap longitudes to [-180, 180] for polygon tests
    lon_wrapped = (((lon2d % 360) + 180) % 360) - 180

    # read Natural Earth land polygons
    shapefile = shpreader.natural_earth(resolution='110m', category='physical', name='land')
    geoms = list(shpreader.Reader(shapefile).geometries())

    # vectorized point-in-polygon
    xv, yv = lon_wrapped.values, lat2d.values
    mask_np = np.zeros(xv.shape, dtype=bool)
    for g in geoms:
        mask_np |= sv.contains(g, xv, yv)

    return xr.DataArray(mask_np, coords={lat_name: lat2d[lat_name], lon_name: lon2d[lon_name]}, dims=lat2d.dims)

def make_land_mask(
    da: xr.DataArray,
    hemisphere: str | None = None,   # "NH", "SH", or None
    resolution: str = "110m",        # "110m", "50m", "10m"
    lon_wrap: bool = True,           # wrap to [-180, 180] for polygon tests
) -> xr.DataArray:
    """
    Return a boolean land mask DataArray with dims/coords aligned to `da`
    and (optionally) restricted to a hemisphere.
    """
    # ---- identify lat/lon coordinate names ----
    lat_name = next(n for n in da.coords if n.lower().startswith("lat"))
    lon_name = next(n for n in da.coords if n.lower().startswith("lon"))

    # ---- build 2D lon/lat on da's grid ----
    lon2d, lat2d = xr.broadcast(da[lon_name], da[lat_name])

    # ---- wrap longitude if desired ----
    if lon_wrap:
        lon2d_test = ((lon2d + 180) % 360) - 180
    else:
        lon2d_test = lon2d

    # ---- hemisphere restriction (applied to mask AFTER land test) ----
    hemi_mask = xr.ones_like(lat2d, dtype=bool)
    if hemisphere is not None:
        hemisphere = hemisphere.upper()
        if hemisphere == "NH":
            hemi_mask = lat2d >= 0
        elif hemisphere == "SH":
            hemi_mask = lat2d < 0
        else:
            raise ValueError("hemisphere must be 'NH', 'SH', or None")

    # ---- Natural Earth land polygons ----
    shapefile = shpreader.natural_earth(
        resolution=resolution, category="physical", name="land"
    )
    geoms = list(shpreader.Reader(shapefile).geometries())

    # ---- compute land mask on numpy arrays ----
    mask_np = np.zeros(lon2d.shape, dtype=bool)
    xv = lon2d_test.values
    yv = lat2d.values
    for g in geoms:
        mask_np |= sv.contains(g, xv, yv)

    mask_da = xr.DataArray(
        mask_np,
        coords={lat_name: da[lat_name], lon_name: da[lon_name]},
        dims=lon2d.dims,  # preserve broadcast dim order
        name="land_mask",
    )

    # Apply hemisphere restriction (preserves dims)
    mask_da = mask_da & hemi_mask

    # ---- broadcast to match da if da has extra dims (time, member, etc.) ----
    mask_da, _ = xr.broadcast(mask_da, da)

    return mask_da