# RECCAP Southern Ocean figures by Cara Nissen

## Figures
-	Figure 2: annual mean $\mathrm{CO}_2$ flux, 1985-2018
-	Figure 3: natural vs. anthropogenic components of $\mathrm{CO}_2$ flux + climate effect
-	Figure 4: winter mean $\mathrm{CO}_2$ flux, 2015-2018
-	Figure 5: summer mean $\mathrm{CO}_2$ flux, 2015-2018
-	Figure 6, part 1: seasonality of $\mathrm{CO}_2$ flux, 1985-2018
-	Figure 6, part 2: maps of season of max. $\mathrm{CO}_2$ uptake, 1985-2018
-	Figure 8: time series 1985-2018 of $\mathrm{CO}_2$ flux
-	Figure 10: trends in $\mathrm{CO}_2$ flux and SST, 1985-2018
-	Figure S1 time series of $\mathrm{CO}_2$ flux in simB, stored trends are loaded into script for Fig. 8 
-	Figure S2: see Fig. 6, part 1, here including MPI model
-	Figure S5: see Figs 2, 4 & 5, zonally averaged $\mathrm{CO}_2$ flux 1985-2018 for individual models
-	Figure S6: see Fig. 3, further decomposition of climate effect 
-	Figure S7: see Fig. 3, further regional decomposition
-	Figure S8: see Fig. 2, here for 2015-2018
-	Figure S9: see Fig. 6, part 2: individual $p\mathrm{CO}_2$-data products
-	Figure S10: see Fig. 6, part 1, further regional decomposition 
-	Figure S12: see Fig. 8, further regional decomposition
-	Figure S14: see Fig. 10, here for 1985-2000
-	Figure S15: see Fig. 10, here for 2001-2018
-	Figure S16: see Fig. 10, trends 1985-2018 in $\mathrm{CO}_2$ flux and SST for individual models
## Data

All scripts require a download of the full RECCAP2 dataset. The full dataset can be downloaded from [https://doi.org/10.5281/zenodo.7990822](https://doi.org/10.5281/zenodo.7990822).

In addition, the notebook [RECCAPv2_SO_Fig10_CO2_flux_trends.ipynb](RECCAPv2_SO_Fig10_CO2_flux_trends.ipynb) requires satellite-based SST data, in particular NOAA’s Extended Reconstructed Sea Surface Temperature, version 5 (ERSSTv5). These data can be downloaded from [https://doi.org/10.7289/V5T72FNM](https://doi.org/10.7289/V5T72FNM) (last access August 30, 2023). 

## Scripts 
-	Figure 2: `RECCAPv2_SO_Fig2_annual_mean_CO2_flux.ipynb`
-	Figure 3: `RECCAPv2_SO_Fig3_bar_plot_anthr_vs_natural.ipynb`
-	Figure 4: `RECCAPv2_SO_Fig4_5_seasonal_mean_CO2_flux_2015_2018.ipynb`
-	Figure 5: `RECCAPv2_SO_Fig4_5_seasonal_mean_CO2_flux_2015_2018.ipynb`
-	Figure 6, part 1: `RECCAPv2_SO_Fig6_seasonal_cycle_part1.ipynb`
-	Figure 6, part 2: `RECCAPv2_SO_Fig6_seasonal_cycle_part2.ipynb`
-	Figure 8: `RECCAPv2_SO_Fig8_time_series.ipynb` (see also script for Fig. S1)
-	Figure 10: `RECCAPv2_SO_Fig10_CO2_flux_trends.ipynb`
-	Figure S1: `RECCAPv2_SO_Suppl_time_series_trends_simB.ipynb`
-	Figure S2: see script for Fig. 6
-	Figure S5: for panel a, see script for Fig. 2; for panels b-c, see script for Figs. 4-5
-	Figure S6: `RECCAPv2_SO_Suppl_bar_plot_climate_effect_only.ipynb`
-	Figure S7: see script for Fig. 3
-	Figure S8: `RECCAPv2_SO_Suppl_annual_mean_CO2_flux_2015_2018_only.ipynb`
-	Figure S9: see script for Fig. 6 (part 2)
-	Figure S10: see script for Fig. 6 (part 1)
-	Figure S12: see script for Fig. 8
-	Figure S14: see script for Fig. 10
-	Figure S15: see script for Fig. 10
-	Figure S16: see script for Fig. 10

To run these scripts, you must have the required packages installed

## Required packages

### Commonly available
-	`os`
-	`numpy`
-	`seawater`
-	`matplotlib`
-	`netCDF4`
-	`datetime`
-	`copy`
-	`cartopy`
-	`scipy`
-	`sklearn`
