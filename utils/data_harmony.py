### some functions
def preprocess(ds):
    """
    Different variable names in time
    """
    if "valid_time" in ds:
        ds = ds.rename({"valid_time": "time"})
    if "valid_time_bnds" in ds:
        ds = ds.rename({"valid_time_bnds": "time_bnds"})
    if "pressure_level" in ds:    
        ds = ds.rename({"pressure_level": "level"})    
    return xr.decode_cf(ds)  # makes time usable if units are valid
