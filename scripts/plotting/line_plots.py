from matplotlib import pyplot as plt


def style_line_subplot(ax, add_zero_line=True, xlim=None):
    import numpy as np
    import pandas as pd
    
    if ax is None:
        ax = plt.gca()
    plt.sca(ax)
    
    plt.xticks(rotation=0, ha='center')
    plt.xlabel('')
    
    if ax.get_ylim()[0] < 0 < ax.get_ylim()[1]:
        ax.axhline(0, color='#CCCCCC', ls='-', lw=0.8, zorder=0)
    
    if xlim is None:
        xlim = ax.get_xlim()
    ax.set_xticks(np.arange('1980', '2020', 5, dtype='datetime64[Y]'))
    ax.set_xticklabels(np.arange(1980, 2020, 5))
    ax.set_xlim(*xlim)
        
    return ax


def plot_ensemble(da, ax=None, dim='variable', draw_stdev=True, color='C0', name='', **kwargs):
    if ax is None:
        ax = plt.gca()
    
    x_dim = list(set(da.dims) - set([dim]))[0]
    x = da[x_dim].values
        
    ens = da
    avg = da.mean(dim)
    
    if draw_stdev:   
        std = da.std(dim)
        upr = avg + std
        lwr = avg - std
        # drawn first to maintain zorder
        ax.fill_between(x, upr, lwr, color=color, alpha=0.25, **kwargs)
        
    ens.plot(c=color, lw=0.5, alpha=0.7, ax=ax, add_legend=False, hue=dim, **kwargs)
    avg.plot(c=color, lw=2.5, alpha=1.0, ax=ax, label=name, **kwargs)
    
    xlim = x.min(), x.max()
    ax = style_line_subplot(ax, add_zero_line=False, xlim=xlim)
    
    return ax