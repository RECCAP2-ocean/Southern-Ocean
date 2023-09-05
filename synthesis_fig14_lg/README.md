# RECCAP-SO: Summary figure

- **Figure 14**: A summary figure of the RECCAP2-Southern Ocean chapter

## Procedure

Map and sections are created in Python as images. The images are then combined in Photoshop. Layers are masked and transformed (warped, flipped, etc.) to create the final image without annotations. Annotations are created in Powerpoint.

For either the Photoshop or Powerpoint files, please contact Luke Gregor (<gregorl@ethz.ch>) directly.

## Data and plotting sections

Data is stored in the `data` folder and contains information for the two sections (Cnet, Cnat), the surface fluxes, and the boundaries of the biomes (smoothed for a nicer looking image). A notebook in the `data` folder creates the plots for the final image.


### Required packages

- `numpy`
- `pandas`
- `xarray`
- `seaborn`
- `cartopy`
- `matplotlib`
- `pooch`