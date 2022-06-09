# Code for Arctic albedo paper

The analysis goes as following

### prepare_MODIS_albedo_and_LSM_python

First, I prepare MCD43A3 broadband albedo with *Albedo_MCD43A3_downoad_mosaic_crop.ipynb* and *Quality_MCD43A2_downoad_mosaic_crop.ipynb* scripts.
The structure of the scripts is the following

1. download MODIS raw albedo HDf files (or QA/snow-flag for each band) from EarthData server via pyModis package for 15th of July each year. For this step, you would need to sign-up to EarthData and insert your log-in info in the first chunk of code.
2. mosaic the HDF files for each year to a single TIFF mosaic and reproject it to the CAVM coordinate system on the 500m resolution
3. crop the resultant mosaic with the CAVM shapefile (above the Arctic treeline)

I also prepare the land surface models (LSMs) NetCDF files: reproject them, convert them to TIFF and crop them via *LSM_NetCDF_to_tiff_reproject_and_crop.ipynb* script. The original LSMs were downloaded for 15th of July for years 2000-2014 of "land-hist" experiment. They are openly available at:
https://esgf-node.llnl.gov/projects/cmip6/

