import xarray as xr
import pandas as pd
import numpy as np
import ocean_data_tools as odt
import pooch


def main():
    dat = (
        xr.open_zarr(
            '../../OceanSODA/ESSD-2020-pH/data-out/'
            'CO2_pCO2m35k68-GBM_TAm30k25-GBM__1985-2018_full_meta.zarr/')
        .rename({'pCO2sea': 'ethGRACER'}))

    xds = downloaders.get_all([dat.ethGRACER])
    regions = get_reccap_regions()

    xds['pco2atm'] = dat.pCO2atm
    xds['sst'] = dat.sst_ostia
    xds['salt'] = dat.salt_soda
    xds['wind_speed'] = dat.windspeed_era5
    xds['press_atm'] = dat.mslp_era5
    xds['ice_frac'] = dat.ice_frac
    xds['regions_reccap2'] = regions.regions_reccap2
    xds['regions_southOcn'] = regions.regions_SO
    xds.pco2sea.attrs = dict(
        standard_name=("surface_partial_pressure_of_"
                       "carbon_dioxide_in_sea_water"),
        units="uatm")

    return xds


class downloaders:
    odt.sparse

    @staticmethod
    def mpiSOMFFN():
        somffn_fname = pooch.retrieve(
            url=("https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/"
                 "0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc"),
            known_hash=None,
            fname='spco2_MPI_SOM-FFN_v2018.nc',
            path='../data-in/interp_products/',
            downloader=pooch.HTTPDownloader(progressbar=True),
        )
        xds = (
            xr.open_dataset(somffn_fname, drop_variables=['date'])
            .prep.rename_to_timelatlon()
            .prep.center_coords_at_0()
            .prep.time_month_day()
        ).spco2_raw
        return xds

    @staticmethod
    def csirML6():
        ml6_fname = pooch.retrieve(
            url=("https://www.nodc.noaa.gov/archive/arc0148/0206205"
                 "/1.1/data/0-data/pCO2_CSIRML6_v2019a.nc"),
            known_hash=None,
            fname='CSIR-ML6_v2019a.nc',
            path='../data-in/interp_products/',
            downloader=pooch.HTTPDownloader(progressbar=True),
        )

        xds = (
            xr.open_dataset(ml6_fname)
            .prep.time_month_day()
        ).pCO2sea_raw
        return xds

    @staticmethod
    def ueaSI():
        fname = pooch.retrieve(
            url=("http://store.pangaea.de/Publications/JonesS-etal_2019/"
                 "statistical_gap-filled_fco2_v1-0b_20180207.nc.zip"),
            known_hash=None,
            fname='UEA-SI_fco2_v2019.nc.zip',
            path='../data-in/interp_products/',
            processor=pooch.Unzip(),
            downloader=pooch.HTTPDownloader(progressbar=True),
        )
        xds = (
            xr.open_mfdataset(fname, combine='nested', concat_dim='time')
            .prep.center_coords_at_0()
            .interp(lat=np.arange(-89.5, 90), lon=np.arange(-179.5, 180))
            .assign_coords(time=pd.date_range('1985', '2018',
                           freq='1MS',
                           closed='left'))
            .roll(lon=20, roll_coords=False)
            .interpolate_na(dim='lon', limit=2)
            .roll(lon=-20, roll_coords=False)
        ).fco2
        return xds

    @staticmethod
    def jenaMLS():
        jenamls_fname = pooch.retrieve(
            url=("http://www.bgc-jena.mpg.de/CarboScope/oc/"
                 "INVERSION/OUTPUT/oc_v1.7_pCO2_daily.nc"),
            known_hash=None,
            fname="JENA_MLS_v17_pCO2_daily.nc",
            path='../data-in/interp_products/',
            downloader=pooch.HTTPDownloader(auth=('CO2inv', '.flux.'),
                                            progressbar=True)
        )

        xds = (
            xr.open_dataset(jenamls_fname, drop_variables=['time'])
            .prep.rename_to_timelatlon()
            .pCO2
            .resample(time='1MS').mean()
            .interp(lat=np.arange(-89.5, 90), lon=np.arange(-179.5, 180))
            .roll(lon=20, roll_coords=False)
            .interpolate_na(dim='lon')
            .roll(lon=-20, roll_coords=False)
        )
        return xds

    @staticmethod
    def lsceFFNNv2(cmems_username=None, cmems_password=None):
        flist = []
        dates = pd.date_range('1985-01', '2019-01', freq='1MS', closed='left')
        for t in dates:
            flist += pooch.retrieve(
                url=(
                    f"ftp://my.cmems-du.eu/Core/MULTIOBS_GLO_BIO_CARBON"
                    f"_SURFACE_REP_015_008/dataset-carbon-rep-monthly/{t:%Y}"
                    f"/dataset-carbon-rep-monthly_{t:%Y%m}15T0000Z_"
                    f"P20191231T1545Z.nc"),
                known_hash=None,
                fname=f"LSCE-FFNNv2-monthly_{t:%Y%m}.nc",
                path='../data-in/interp_products/LSCE-FFNNv2-monthly',
                downloader=pooch.FTPDownloader(
                    username=cmems_username,
                    password=cmems_password,
                    progressbar=True)
            ),

        xds = xr.open_mfdataset(flist, combine='nested', concat_dim='time')
        xds = (xds
               .prep.rename_to_timelatlon()
               .prep.center_coords_at_0()
               .prep.time_month_day()
               .spco2 * 9.86923)

        return xds

    @staticmethod
    def jmaMLR():
        flist = []
        dates = pd.date_range('1990-01', '2019-01', freq='1AS', closed='left')
        for t in dates:
            flist += pooch.retrieve(
                url=(f"http://www.data.jma.go.jp/gmd/kaiyou/data/"
                     f"english/co2_flux/grid/JMA_co2map_{t:%Y}.ZIP"),
                known_hash=None,
                fname=f"JMA-MLR_co2map_{t:%Y}.nc",
                path='../data-in/interp_products/JMA-MLR-monthly',
                processor=pooch.Unzip(),
                downloader=pooch.HTTPDownloader(progressbar=True))

        xds = (
            xr.open_mfdataset(flist,
                              decode_times=False,
                              concat_dim='time',
                              combine='nested')
            .prep.rename_to_timelatlon()
            .prep.center_coords_at_0()
            .assign_coords(time=pd.date_range('1990-01', '2019-01',
                                              freq='1MS',
                                              closed='left'))
        ).pCO2s
        return xds

    @staticmethod
    def get_funcs():
        obj = downloaders
        funcs = [v for k, v in obj.__dict__.items()
                 if callable(getattr(obj, k))]

        funcs = [f for f in funcs if not f.__name__.startswith('get')]
        return funcs

    @staticmethod
    def get_all(xds=[]):

        funcs = downloaders.get_funcs()
        names = []
        if len(xds) > 0:
            for item in xds:
                names += item.name,
        for func in funcs:
            print(func.__name__)
            names += func.__name__,
            xds += func().load(),
        xds = (
            xr.concat(xds, dim='meth')
            .assign_coords(meth=names)
            .to_dataset(name='pco2sea'))
        xds['mask'] = xds.pco2sea.notnull().all('meth')

        return xds


def get_reccap_regions():
    so_regions = pooch.retrieve(
        url=(
            "https://github.com/RECCAP2-ocean/shared-resources/raw/"
            "master/regions/reccap2ocean_SOsubregions.nc"),
        known_hash=None,
        fname='reccap2_regions_SOsub.nc',
        path='../data-in/regions/'
    )

    reccap2_regions = pooch.retrieve(
        url=(
            "https://github.com/RECCAP2-ocean/shared-resources/raw/"
            "master/regions/reccap2ocean_regions.nc"),
        known_hash=None,
        fname='reccap2_regions.nc',
        path='../data-in/regions/',
    )

    regions = (
        xr.open_dataset(reccap2_regions, drop_variables=['option1', 'option2'])
        .rename({
            'option3': 'regions_reccap2',
            'fay_mckinley': 'regions_fayMckinley',
            'woa_regions': 'regions_woa'})
    )

    regions['regions_SO'] = xr.open_dataset(so_regions).so_subregions

    return regions


if __name__ == '__main__':
    pass
