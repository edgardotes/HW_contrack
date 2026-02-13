## Climate indices
import numpy as np
import xarray as xr

def count_frost_days(tn, threshold=0.0):
    '''
    Annual count of days when TN (daily minimum temperature) < 0oC.

    Let TNijbe daily minimum temperature on day i in year j. Count the number of days where:

    TNij < 0oC.
    '''
    # Clean data (handle any non-finite values safely)
    tn_clean = tn.where(np.isfinite(tn))

    # Frost day mask and count
    frost_days = (tn_clean < threshold)
#    return frost_days.sum(dim='time', skipna=True)
    # Convert boolean to float (True=1.0, False=0.0), and sum over time
    return frost_days.astype(float).sum(dim='time', skipna=True)

def count_summer_days(tm, threshold=25.0):
    '''
    SU, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25oC.

    Let TXij be daily maximum temperature on day i in year j. Count the number of days where:
    TXij > 25oC.
    '''
    # Clean data (handle any non-finite values safely)
    tm_clean = tm.where(np.isfinite(tm))

    # Assume °C by default, adjust only if unit says Kelvin
 #   if tn_clean.attrs.get('units', '').lower() in ['k', 'kelvin']:
 #       threshold += 273.15

    # Frost day mask and count
    summer_days = (tm_clean > threshold)
    return summer_days.sum(dim='time', skipna=True)

def count_icing_days(tm, threshold=0.0):
    '''
    ID, Number of icing days: Annual count of days when TX (daily maximum temperature) < 0oC.

    Let TXijbe daily maximum temperature on day i in year j. Count the number of days where:

    TXij < 0oC.
    '''
    # Clean data (handle any non-finite values safely)
    tm_clean = tm.where(np.isfinite(tm))

    # Assume °C by default, adjust only if unit says Kelvin
 #   if tn_clean.attrs.get('units', '').lower() in ['k', 'kelvin']:
 #       threshold += 273.15

    # Frost day mask and count
    icing_days = (tm_clean < threshold)
    return icing_days.sum(dim='time', skipna=True)
#def compute_TXX(tmax):
#    return tmax.groupby('time.year').max(dim='time')

def compute_dtr(tx, tn):
    """
    Compute Daily Temperature Range (DTR) as the monthly mean of (TX - TN).

    Parameters:
        tx (xarray.DataArray): Daily maximum temperature (TX), with 'time' dimension.
        tn (xarray.DataArray): Daily minimum temperature (TN), same dimensions as tx.

    Returns:
        xarray.DataArray: Monthly mean daily temperature range.
    """
    # Ensure both inputs align in time and space
    tx, tn = xr.align(tx, tn)

    # Compute daily range
    dtr_daily = tx - tn

    # Compute monthly mean DTR
    dtr_monthly = dtr_daily.resample(time='ME').mean(dim='time', skipna=True)

    # Add metadata
    dtr_monthly.name = 'DTR'
    dtr_monthly.attrs.update({
        'long_name': 'Monthly mean daily temperature range',
        'units': tx.attrs.get('units', '°C'),
        'comment': 'Mean of daily (TX - TN), where TX is max and TN is min temperature.',
        'cell_methods': 'time: mean (of daily max - min)',
    })

    return dtr_monthly
def get_max_dry_spell(seq):
    """
    Get the length of the longest consecutive run of True (dry days) in a 1D boolean array.
    """
    return max((sum(1 for _ in group) for val, group in groupby(seq) if val), default=0)


def compute_CDD(precip, threshold=1.0):
    """
    Compute the maximum length of dry spells (CDD) using groupby logic.
    
    Parameters:
        precip (xarray.DataArray): Daily precipitation with time dim.
        threshold (float): Dry day threshold in mm (default: 1.0).
    
    Returns:
        xarray.DataArray: CDD values (max dry spell length) at each grid point.
    """
    # Get dry day mask
    dry_days = (precip < threshold)

#    dry_days=dry_days.astype(float)

    # Apply the get_max_dry_spell function along 'time' using xarray.apply_ufunc
    cdd = xr.apply_ufunc(
        get_max_dry_spell,
        dry_days,
        input_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[int]
    )

    cdd.name = 'CDD'
    cdd.attrs.update({
        'long_name': 'Maximum number of consecutive dry days',
        'units': 'days',
        'comment': 'CDD is the length of the longest dry spell (RR < 1 mm)',
    })

    return cdd,dry_days

def compute_CWD(precip, threshold=1.0):
    """
    Compute the maximum length of wet spells (CWD) using groupby logic.
    
    Parameters:
        precip (xarray.DataArray): Daily precipitation with time dim.
        threshold (float): Dry day threshold in mm (default: 1.0).
    
    Returns:
        xarray.DataArray: CWD values (max dry spell length) at each grid point.
    """
    # Get dry day mask
    dry_days = (precip > threshold)

#    dry_days=dry_days.astype(float)

    # Apply the get_max_dry_spell function along 'time' using xarray.apply_ufunc
    cwd = xr.apply_ufunc(
        get_max_dry_spell,
        dry_days,
        input_core_dims=[['time']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[int]
    )

    cwd.name = 'CWD'
    cwd.attrs.update({
        'long_name': 'Maximum number of consecutive wet days',
        'units': 'days',
        'comment': 'CWD is the length of the longest wet spell (RR > 1 mm)',
    })

    return cwd,wet_days

def compute_precip_exceedance(precip, start='2010', end='2023', bin_width=0.01, max_precip=1000):
    """
    Compute exceedance frequency of daily precipitation.
    
    Parameters:
        precip (xarray.DataArray): Daily precipitation (mm/day), must have 'time' dim.
        start (str): Start year.
        end (str): End year.
        bin_width (float): Width of precipitation bins (mm/day).
        max_precip (float): Upper limit of precipitation to consider (mm/day).
    
    Returns:
        bins (np.ndarray): Precipitation bin edges.
        exceedance_freq (np.ndarray): Exceedance frequency (%) for each bin.
    """
    # Select time slice
    p = precip.sel(time=slice(start, end))
    
    # Flatten across time and space
    values = p.values.flatten()
    values = values[np.isfinite(values)]  # Remove NaNs

    # Define bins
    bins = np.arange(0, max_precip + bin_width, bin_width)

    # Histogram (counts of values in each bin)
    hist, bin_edges = np.histogram(values, bins=bins)

    # Cumulative sum (reverse): exceedance count
    exceedance_counts = np.cumsum(hist[::-1])[::-1]

    # Convert to percentage
    exceedance_freq = 100 * exceedance_counts / exceedance_counts[0]

    return bins[:-1], exceedance_freq

def get_trend(xx: np.ndarray, yy: np.ndarray) -> (xr.DataArray, xr.DataArray, xr.DataArray):
    """Fit linear trend and p-value.

    Parameters
    ----------
    xx : np.ndarray, shape (N,)
    yy : np.ndarray, shape (N,)

    Returns
    -------
    tuple of DataArray : (slope, intercept, pvalue)
    """
    results = linregress(xx, yy, alternative='greater')  # only look for statistically significant positive trends
    return (
        xr.DataArray(results.slope, name='slope'),
        xr.DataArray(results.intercept, name='intercept'),
        xr.DataArray(results.pvalue, name='pvalue'),
    )
