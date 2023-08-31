def plot_map(ds = None, dsg = None, var = 'SiO3', fig = None, pos = 111, contour = False, cont_levs = None,
            colbar_ori = "vertical", colbar_loc = [0.95, 0.1, 0.02, 0.7], cbar_ticks = None,
            LAT_S = -90, LAT_N = -30, LON_W = -180, LON_E = 180, proj = 'SouthPolarStereo', proj_shape = "round",
            colbar = True, lev = -1, timestep = 0, varc = None, clevs = None, extend = "neither", ncols_per_int = 2,
            hatch_mask = None, hatches = '////', grid = 'roms', cmap_sym = False, title = None, clabel = True,
            scale = 1.0, vmin = None, vmax = None, ncols = 10, clrmap = 'viridis', cmap_lib = 'cmocean',
            cbar_label = 'Pg C/y', lonname = 'lon', latname = 'lat', grid_label_lon = None, grid_roms = "rho", axes = None,
            label_basins = False, basin_lons = [24, 145, -70], label_lines = True, label_grid_lines = True):

    """Plots 2d data on a cartopy map, with optional hatching and/or contours overlain. 
    
    Keyword arguments:
    Required:
    ds: An xarray dataset containing the variable to be plotted.
    dsg: An xarray dataset containing the grid on which the variable is stored (can be the same dataset as above).
    fig: Figure instance for plotting
    var (str): String identifying the variable to plot in pcolormesh
    pos (int):  a three digit integer, where the first digit is the number of rows, the second the number of columns, and the third the index of the subplot (default: 111).
    axes (axis instance): Provide the plotting axis directly, *instead* of the pos argument. 

    Optional:
    grid (str): {"roms", "standard"} (default: "roms"). Grid of the input data. 
    proj (str): Cartopy map projection class. Must match exactly the class name, e.g. 'Mercator', 'SouthPolarStereo'
    lev (int): Depth level index for 4D data (default: -1)
    timestep (int): Time index for 3D/4D data (default: 0)

    varc (str): Name of variable to be used for overlain contours, must be part of Dataset "ds". Optional.
    clevs (array-like): List or array of contour levels.

    hatch_mask (array-like): Array used to overlay hatching. Must be equal to 1 where hatching is desired. 
    hatches (str): Hatching style according to matplotlib. Default is '////'.

    colbar (bool): Draw a colorbar or not (default: True)
    colbar_loc (list): Colorbar axis. As for matplotlib.pyplot.figure.add_axes 
    colbar_ori (str): {"horizontal", "vertical"}. As for matplotlib.pyplot.figure.colorbar 
    cmap_lib (str): {"cmocean", "matplotlib"} Color map library to use.
    clrmap (str or cmap): Name of matplotlib or cmocean colormap for pcolormesh plot (default cmocean: cmo.haline; default matplotlib: 'viridis')
    ncols (int): Number of colors in clrmap (default: 10)
    cmap_sym (bool): Whether the colormap should be symmetrical around zero or not (default: False)
    vmin, vmax (scalar): As for matplotlib.pyplot.pcolormesh
    cbar_label (str): Label for colorbar.
    extend (str): {"both", "max", "min", "neither"} (default: "neither"). As for matplotlib.pyplot.colorbar

    scale (scalar): Constant which "var" is multiplied by before plotting (default: 1.0)
    lonname (str): Name of longitude variable if "grid" is "standard" (default: 'lon')
    latname (str): Name of latitude variable if "grid" is "standard" (default: 'lat')
    grid_label_lon (float): Longitude to label 60S and 40S if proj = 'SouthPolarStereo'
    """

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm
    import numpy as np

    import xarray as xr
    import pandas as pd
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.path as mpath
    import sys
    import cmocean
    from matplotlib.colors import LinearSegmentedColormap


    def fmt(x):
        """ Removes trailing zeros for contour labelling"""
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return s

    proj_obj = eval('ccrs.'+proj+'()')

    #fig = plt.figure(figsize=[10, 10])
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16)
    if axes is None:
        ax = fig.add_subplot(pos, projection= proj_obj)
    else:
        ax = fig.add_subplot(axes, projection= proj_obj)
    cbar_loc = colbar_loc
    #
    #cmap.set_bad('grey')
    if cmap_lib == 'matplotlib':
        cmap = cm.get_cmap(clrmap, ncols)
    elif cmap_lib == 'cmocean':
        if clrmap == 'viridis':
            clrmap = cmocean.cm.haline
        cdict = cmocean.tools.get_dict(clrmap, N=ncols)
        cmap = LinearSegmentedColormap('Cmap', segmentdata=cdict, N=ncols)
        if '_r' in clrmap.name:
            cmap = cmap.reversed()

    if grid == 'roms':
        if grid_roms == "rho":
            lonn = np.where(dsg.lon_rho > 360, dsg.lon_rho - 360, dsg.lon_rho)
            lat = dsg.lat_rho.values
            lonn = lonn[:,1:-1]
            lat = lat[:,1:-1]
            lon = np.roll(lonn, 24*4, axis = 1)
        elif grid_roms == "u":
            lonn = np.where(dsg.lon_u > 360, dsg.lon_u - 360, dsg.lon_u)
            lat = dsg.lat_u.values
            lonn = lonn[:,0:-1]
            lat = lat[:,0:-1]
            lon = np.roll(lonn, (24*4) - 1, axis = 1)
        
         
        
    elif grid == 'standard':
        lon = dsg[lonname]
        lat = dsg[latname]
        if lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
    else:
        print('ERROR: Grid not yet supported! Pass either roms or standard lon-lat grid.')
    
    ndims = len(ds[var].shape)
    if ndims == 4:
        varr = ds[var][timestep,lev,:,:]*scale
    elif ndims == 3:
        varr = ds[var][timestep,:,:]*scale
    elif ndims == 2:
        varr = ds[var]*scale
    
    if varc is not None:
        varc = ds[varc]
    # remove east-west ghost points
    if grid == 'roms':
        if grid_roms == "rho":
            varr = varr[:,1:-1]
        elif grid_roms == "u":
            varr = varr[:,0:-1]
        if varc is not None:
            if grid_roms == "rho":
                varc = varc[:,1:-1]
            elif grid_roms == "u":
                varc = varc[:,0:-1]
        
        if proj == 'SouthPolarStereo':
            varr = np.roll(varr, 24*4 +1, axis = 1)
            if varc is not None:
                varc = np.roll(varc, 24*4 +1, axis = 1)
            if hatch_mask is not None:
                if grid_roms == "rho":
                    hatch_mask = hatch_mask[:,1:-1]
                elif grid_roms == "u":
                    hatch_mask = hatch_mask[:,0:-1]
                hatch_mask = np.roll(hatch_mask, 24*4, axis = 1)

    # compute default vmin, vmax
    if vmin is None:
        vmin = np.nanpercentile(varr, 1)
    if vmax is None:
        vmax = np.nanpercentile(varr, 99)

    if cmap_sym == True:
        if np.abs(vmin) > vmax:
            vmax = np.abs(vmin)
        elif np.abs(vmax) > np.abs(vmin):
            vmin = -1*vmax 

    # make values outside vmin, vmax be the last/first color of the cmap
    # adapted from https://stackoverflow.com/questions/48613920/use-of-extend-in-a-pcolormesh-plot-with-discrete-colorbar
    if (extend == "both") or (extend == "max") or (extend == "min"):
        # create boundaries between vmin and vmax, to have ncols intervals
        boundaries = np.linspace(vmin, vmax, np.int((ncols / ncols_per_int)+1) )
        # create list of colors from colormap
        if extend == "max":
            if cmap_lib == 'matplotlib':
                cmap_tmp = plt.cm.get_cmap(clrmap, len(boundaries))
            elif cmap_lib == 'cmocean':
                cdict = cmocean.tools.get_dict(clrmap, N = len(boundaries))
                cmap_tmp = LinearSegmentedColormap('Cmap', segmentdata=cdict, N=len(boundaries))
                if '_r' in clrmap.name:
                    cmap_tmp = cmap.reversed()
            colors = list(cmap_tmp(np.arange(len(boundaries))))
            cmap = mpl.colors.ListedColormap(colors[:-1], "")
            cmap.set_over(colors[-1])
        elif extend == "min":
            if cmap_lib == 'matplotlib':
                cmap_tmp = plt.cm.get_cmap(clrmap, len(boundaries))
            elif cmap_lib == 'cmocean':
                cdict = cmocean.tools.get_dict(clrmap, N = len(boundaries))
                cmap_tmp = LinearSegmentedColormap('Cmap', segmentdata=cdict, N=len(boundaries))
                if '_r' in clrmap.name:
                    cmap_tmp = cmap.reversed()
            colors = list(cmap_tmp(np.arange(len(boundaries))))
            cmap = mpl.colors.ListedColormap(colors[1:], "")
            cmap.set_under(colors[0])
        elif extend == "both":
            if cmap_lib == 'matplotlib':
                cmap_tmp = plt.cm.get_cmap(clrmap, ncols+2)
            elif cmap_lib == 'cmocean':
                cdict = cmocean.tools.get_dict(clrmap, N = ncols+2)
                cmap_tmp = LinearSegmentedColormap('Cmap', segmentdata=cdict, N=ncols+2)
                if '_r' in clrmap.name:
                    cmap_tmp = cmap_tmp.reversed()
            colors = list(cmap_tmp(np.arange(ncols+2)))
            cmap = mpl.colors.ListedColormap(colors[1:-1], "")
            cmap.set_over(colors[-1])
            cmap.set_under(colors[0])

    
    # Limit the map to -30 degrees latitude and below.
    ax.set_extent([LON_W, LON_E, LAT_S, LAT_N], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor = 'grey', zorder = 10)


    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    if (proj == 'SouthPolarStereo') and (proj_shape == "round"):
        theta = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)
    elif (proj == 'SouthPolarStereo') and (proj_shape == "square"):
        pass

    if label_grid_lines:
        gl = ax.gridlines(crs = ccrs.PlateCarree())

    if contour:
        pcol = ax.contourf(lon, lat, varr, transform = ccrs.PlateCarree(), cmap = cmap,  vmin = vmin, vmax = vmax,
                         levels = cont_levs, zorder = -1, extend = extend)
    else:
        pcol = ax.pcolormesh(lon, lat, varr, transform = ccrs.PlateCarree(), cmap = cmap, 
                         vmin = vmin, vmax = vmax, zorder = -1, shading='auto')

    if clevs is not None:
        cont = ax.contour(lon, lat, varc, levels = clevs, transform = ccrs.PlateCarree(), colors = 'k', zorder = 3 )
        if clabel:
            ax.clabel(cont, clevs, inline = True, fontsize = 10, fmt = fmt)

    if hatch_mask is not None:
        hatch = ax.fill_between([np.min(lon),np.max(lon)], np.min(lat), np.max(lat), transform = ccrs.PlateCarree(),
                                hatch=hatches, color="none", edgecolor='black', zorder = 1)

        ax.pcolormesh(lon, lat, np.where(hatch_mask == 1, varr, np.nan), zorder = 2, cmap = cmap,  shading='auto',
                 vmin = vmin, vmax = vmax, transform = ccrs.PlateCarree())

    if colbar == True:
        cb_ax = fig.add_axes(cbar_loc)
        cbar = fig.colorbar(pcol, cax=cb_ax, extend = extend, orientation = colbar_ori)
        if cbar_ticks:
            cbar.set_ticks(cbar_ticks)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label(cbar_label, fontsize = 16)

    if title is not None:
        ax.set_title(title, fontsize = 16, fontweight = 'bold')

    if label_basins:
        # plt_crr = ccrs.PlateCarree()
        ax.plot([basin_lons[0],basin_lons[0]],[-80,-30], color = 'k', lw = 2, transform = ccrs.PlateCarree())
        ax.plot([basin_lons[1],basin_lons[1]],[-80,-30], color = 'k', lw = 2, transform = ccrs.PlateCarree())
        ax.plot([basin_lons[2],basin_lons[2]],[-70,-30], color = 'k', lw = 2, transform = ccrs.PlateCarree())



    # manually label the grid lines
    if (proj == 'SouthPolarStereo') and (label_lines == True):
        ax.text(0, -34.5, '0$^\degree$', transform = ccrs.PlateCarree(), fontsize = 16)
        ax.text(60, -29, '60$^\degree$E', transform = ccrs.PlateCarree(), fontsize = 16)
        ax.text(120, -29, '120$^\degree$E', transform = ccrs.PlateCarree(), fontsize = 16)

        ax.text(180, -31, '180$^\degree$', transform = ccrs.PlateCarree(), fontsize = 16)
        ax.text(-60, -25, '60$^\degree$W', transform = ccrs.PlateCarree(), fontsize = 16)
        ax.text(-120, -24, '120$^\degree$W', transform = ccrs.PlateCarree(), fontsize = 16)

        if grid_label_lon is None:
            ax.text(0, -41, '40$^\degree$S', transform = ccrs.PlateCarree(), fontsize = 16)
            ax.text(0, -61, '60$^\degree$S', transform = ccrs.PlateCarree(), fontsize = 16)
        else:
            ax.text(grid_label_lon, -41, '40$^\degree$S', transform = ccrs.PlateCarree(), fontsize = 12)
            ax.text(grid_label_lon, -61, '60$^\degree$S', transform = ccrs.PlateCarree(), fontsize = 12)

