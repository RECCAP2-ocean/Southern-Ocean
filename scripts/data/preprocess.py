"""
Intended for data to be used with:
    `xr.open_mfdataset(fname, preprocess=reccap2_dataprep)`

Ensures that data is formatted so that it is formatted to match XYZT
    lon (0.5 : 360)
    lat (-89.5 : 90)
    depth (positive downward)
    time (days since 1980-01-01)
"""

import xarray as xr
import pandas as pd
import numpy as np


def add_history(xds, message):
    time = pd.Timestamp.today().strftime('%Y-%m-%d')
    message = message.replace("'", "").replace('"', "")
    
    hist = xds.attrs.get('history', '').strip()
    hist = hist if hist.endswith(';') or hist == '' else hist + '; '
    hist += f' [RECCAP @ {time}] {message}; '
    
    xds = xds.assign_attrs(history=hist.strip())
    
    return xds


def try_run(f):
    import warnings
    import functools
    
    @functools.wraps(f)
    def wrapper(xds):
        try:
            return f(xds)
        except Exception as e:
            warnings.warn(f'{f.__name__} could not be run: {str(e)}')
            return xds
    return wrapper


@try_run
def rename_coords(
    xds,
    lat=['latitude', 'Lat', 'ylat'],
    lon=['longitude', 'Lon', 'xlon'],
    time=['Time', 'tmnth'],
    region=['reg', 'regions'],
    **kwargs
):
    
    coord_names = dict(
        lat=lat, 
        lon=lon,
        time=time,
        **kwargs,
    )
    
    rename = {}
    for key in xds.coords:
        for new, old in coord_names.items():
            if key in old:
                rename[key] = new
    
    xds = xds.rename(**rename)
    
    if rename != {}:
        xds = add_history(xds, f'renamed variables to match protocol {str(rename)}')
    
    return xds


@try_run
def order_lon_360(xds):
    """assumes lon is in xds"""
    if 'lon' not in xds:
        return xds
    
    lon = xds.lon.copy().values
    xds = xds.assign_coords(lon=lambda a: a.lon % 360).sortby('lon')
    
    if any(lon != xds.lon.values):
        xds = add_history(xds, 'Flipped longitudes 0:360')

    return xds 


@try_run
def transpose(xds):
    reccap_order = ["lon", "lat", "depth", "time"]
    coords = list(xds.dims)
    
    order = []
    for key in reccap_order:
        if key in coords:
            order += key,
            coords.remove(key)
    order += coords
    xds = xds.transpose(*order)
    xds = add_history(xds, f'Transposed dimensions to ({str(order)[1:-1]})')

    return xds


@try_run
def time_decoder(xds):
    if 'time' not in xds:
        return xds

    t = xds.time.values.astype(int)

    # nies workaround (they use seconds since 1980)
    if t[0] > 10000:
        t = (t / 86400).astype(int)
    
    xds['time'] = xr.DataArray(
        data=t, dims=('time',), 
        coords={'time':t}, 
        attrs=dict(units='days since 1980', calendar='gregorian'))
    xds = xr.decode_cf(xds)
    
    return xds


@try_run
def center_time_on_15th(xds):
    if 'time' not in xds:
        return xds
    
    dt = np.timedelta64(14, 'D')
    t = xds.time.values.astype('datetime64[M]') + dt
    
    xds['time'] = xr.DataArray(
        data=t, dims=['time'], coords={'time': t})
    
    xds = xds.sel(time=slice('1980', '2018'))
    xds = add_history(xds, f'Time set to 15th of each month and 198[0/5]-2018')
    
    
    return xds
    

def reccap2_dataprep(decode_times=False):
    def reccap2_dataprep(ds):
        ds = rename_coords(ds)
        if decode_times:
            ds = time_decoder(ds)
        ds = center_time_on_15th(ds)
        ds = order_lon_360(ds)
        ds = transpose(ds)
        return ds
    return reccap2_dataprep

