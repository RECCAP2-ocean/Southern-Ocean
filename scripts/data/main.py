from .utils import read_reccap2_products as _read_reccap2_products
from .download import retrieve_files as download, flatten_list as _flatten, read_catalog
from .preprocess import reccap2_dataprep 

data_surface_CO2_vars = [    
    { # variables
        'fgco2_reg': "fgco2_reg",
        'fgco2_glob': "fgco2_glob",
        'fgco2': '(fgco2)\w(?!glob|reg)',
        'spco2': "spco2",
        'pco2atm': 'p.?co2.?atm',
        'Kw': "Kw",
        'alpha': "(alpha)\w(?!.*skin)",
#         'alpha_skin': "alpha_skin",
#         'alpha_subskin': "alpha_subskin",
        'fice': 'siconc|fice',
#         'talk': "talk",
#         'ph': "(?!al)(ph)(?!a)",
#         'tos': "tos",
#         'sos': "sos",
#         'dissic': "dissic",
        'area': "area"},
    { # data products
        'UOEX': 'UOEX_Wat20',
        'NIES': 'NIES',
        'ETHZ': 'OceanSODAETHZ',
        'CSIR': 'CSIRML6',
        'CMEMS': 'CMEMS-LSCE-FFNN'}
]


def reccap_surface_CO2(
    data_dict, 
    data_surface_CO2_vars=data_surface_CO2_vars, 
    verbose=True, 
    **kwargs
):
    """
    kwargs are passed to data.downloads.retrieve_files
    """
    flist = download(**data_dict, verbose=verbose, **kwargs)
    flist = _flatten(flist)
    obj = _read_reccap2_products(flist, data_surface_CO2_vars)
    
    return obj


def soccom_float_liar(data_dict, verbose=True, **kwargs):
    from .non_reccap_data import grid_soccom_argo_float
    
    flist = _flatten(download(**data_dict, verbose=verbose, **kwargs))
    flist = [f for f in flist if f.endswith('.nc')]
    
    variables = [
        'Temperature', 'Salinity', 
        'pCO2_LIAR', 'TALK_LIAR', 'DIC_LIAR', 'TALK_LIAR_QFA', 'DIC_LIAR_QFA']

    liar = grid_soccom_argo_float(flist, keep_vars=variables, agg_funcs=['mean'])

    qc = (liar.DIC_LIAR_QFA == 0) & (liar.TALK_LIAR_QFA == 0)
    liar = liar.where(qc).drop(['DIC_LIAR_QFA', 'TALK_LIAR_QFA'])

    liar = reccap2_dataprep()(liar)
    
    liar = liar.rename(
        Salinity='so',
        Temperature='thetao',
        pCO2_LIAR='pco2',
        TALK_LIAR='talk',
        DIC_LIAR='dissic',
    )

    return liar


def socat2020(
    data_dict, 
    variables=['fco2_ave_unwtd', 'fco2_std_unwtd'], 
    verbose=True, 
    **kwargs
):
    """kwargs are passed to scripts.data.download.retrieve_files"""
    import xarray as xr
    socat = xr.open_mfdataset(
        download(**data_dict, verbose=verbose, **kwargs), 
        preprocess=reccap2_dataprep(decode_times=False)
    )[variables]
    
    return socat


def southern_ocean_mask():
    from scripts.data.masks import get_southern_ocean_subregions
    return get_southern_ocean_subregions()[['subregions', 'names']]