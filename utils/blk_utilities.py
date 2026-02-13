#Blocking utilities
import numpy as np

def reduce_to_1D(ds, latitude_range, latitude_name='Latitude',
                 time_mean=True, time_name='Time'):
    """
    TODO
    """
    if isinstance(latitude_range, (int, float)):
        latitude_range = [latitude_range, latitude_range]
    elif len(latitude_range) == 1:
        latitude_range = [latitude_range[0], latitude_range[0]]
    elif len(latitude_range) != 2:
        errmsg = ' '.join(['latitude_range has to be float or list of',
                           'one or two floats and not {}']).format(
                               latitude_range)
        raise ValueError(errmsg)
    lats = ds[latitude_name].data
    lats = lats[(lats >= latitude_range[0]) & (lats <= latitude_range[1])]
    ds = ds.sel(**{latitude_name: lats})
    ds = (ds.sum(latitude_name) > 0).astype(int)

    if time_mean:
        ds = ds.mean(time_name)
    return ds
    
#def reduce_to_1D(self, latitude_range, time_mean=True, inplace=True):
#        ds = ut.reduce_to_1D(self.ds,
#                             latitude_range=latitude_range,
#                             latitude_name=self._latitude_name,
#                             time_mean=time_mean,
#                             time_name=self._time_name)
#        if not inplace:
#            return ds
#        self.ds = ds
#        
#blk.reduce_to_1D([0, 75])

def calculate_smoothed_field(
    data,
    passes,
    weights=np.array([[0, 1, 0], [1, 2, 1], [0, 1, 0]]),
    mode="wrap",
    *args,
    **kwargs
):
    """
    Calculate smoothed field based on a two-dimensional weight kernel
    and multiple smoothing passes. Default weight kernel is a 3x3
    5-point smoothing with double-weighted centre. The arguments
    "weight" and "mode" must be accepted by scipy.ndimage.convolve.
    Values at the latitude border are always set to NaN.
    Dimension names ("time_name", "lon_name", "lat_name"), size ("ntime", "nlon", "nlat")
    and resolution ("dlon", "dlat") can be passed as key=value arguments.

    Parameters
    ----------
        data : xarray.DataArray
            data to smooth
        passes : int or float
            number of smoothing passes of the 5-point smoothing
        weigths : array_like, optional
            array of weight, two-dimensional
            (see scipy.ndimage.convolve function)
        mode : string, optional
            defines how the array is extended at boundaries
            (see scipy.ndimage.convolve function)

    Returns
    -------
        smoothed data: xarray.DataArray
            Data containing the smoothed field
    """

    # perform smoothing
    smoothed = []
    for step in data[kwargs["time_name"]]:
        temp = data.sel({kwargs["time_name"]: step})
        for p in range(passes):
            temp = ndimage.convolve(temp, weights=weights, mode=mode) / np.sum(weights)

        # set latitude border values to nan
        border_size = int(weights.shape[0] / 2 + 0.5)
        temp[np.arange(-border_size, border_size), :] = np.nan

        smoothed.append(temp)

    # define DataArray
    da = xr.DataArray(
        smoothed,
        coords=[
            data[kwargs["time_name"]],
            data[kwargs["lat_name"]],
            data[kwargs["lon_name"]],
        ],
    )

    # set name
    da.name = "smooth_" + data.name

    # assign attributes
    da = da.assign_attrs(data.attrs)
    da.attrs["smooth_passes"] = passes

    return da

