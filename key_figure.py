# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: gmt
#     language: python
#     name: gmt
# ---

# %% [markdown]
# # **ICESat-2 detected active subglacial lakes**
#
# Making pretty maps of active subglacial lakes in Antarctica,
# a companion jupyter notebook to https://github.com/weiji14/deepicedrain.
# Uses [PyGMT](https://www.pygmt.org) for illustration,
# heavily inspired by https://github.com/mrsiegfried/Venturelli2020-GRL

# %%
import os
from glob import glob

import geopandas as gpd
import numpy as np
import pandas as pd
import pooch
import pygmt

# %% [markdown]
# # Get data files
#
# Data for background basemaps:
#
# - Haran, T. M., Bohlander, J., Scambos, T. A., Painter, T. H., & Fahnestock, M. A. (2014). MODIS Mosaic of Antarctica 2008-2009 (MOA2009) Image Map. U.S. Antarctic Program Data Center (USAP-DC), via National Snow and Ice Data Center (NSIDC). https://doi.org/10.7265/N5KP8037
# - Depoorter, M. A., Bamber, J. L., Griggs, J. A., Lenaerts, J. T. M., Ligtenberg, S. R. M., van den Broeke, M. R., & Moholdt, G. (2013). Antarctic masks (ice-shelves, ice-sheet, and islands), link to shape file (p. 15.8 MBytes) [Application/zip]. PANGAEA - Data Publisher for Earth & Environmental Science. https://doi.org/10.1594/PANGAEA.819147
# - Mouginot, J., Rignot, E., & Scheuchl, B. (2019). MEaSUREs Phase Map of Antarctic Ice Velocity, Version 1 [Data set]. NASA National Snow and Ice Data Center DAAC. https://doi.org/10.5067/PZ3NJ5RXRH10

# %%
# Using Quantarctica3 from https://www.npolar.no/quantarctica/
datafold: str = os.getenv("DATAHOME") or os.path.abspath("Quantarctica3")
os.makedirs(name=datafold, exist_ok=True)

# %%
# MODIS Mosaic of Antarctica
moa_no_nan: str = pooch.retrieve(
    url="ftp://ftp.nsidc.org/pub/DATASETS/nsidc0593_moa2009/geotiff/moa750_2009_hp1_v01.1.tif.gz",
    known_hash="90d1718ea0971795ec102482c47f308ba08ba2b88383facb9fe210877e80282c",
    path=f"{datafold}/SatelliteImagery/MODIS",
    processor=pooch.Decompress(name="moa750_2009_hp1_v1.1.tif"),
)
moa = f"{datafold}/SatelliteImagery/MODIS/moa750_2009_hp1_v01.1.tif"
try:
    assert os.path.exists(path=moa)
except AssertionError:
    with pygmt.clib.Session() as lib:
        # !gmt grdmath $moa_no_nan 0 NAN = $moa
        lib.call_module(module="grdmath", args=f"{moa_no_nan} 0 NAN = {moa}")

# %%
# Scripps Grounding Line
shapefiles: list = pooch.retrieve(
    url="https://epic.awi.de/id/eprint/33781/1/Antarctica_masks.zip",
    known_hash="e4c5918240e334680aed1329f109527efd8f43b6a277bb8e77b66f84f8c16619",
    fname="groundingline",
    path=f"{datafold}/Miscellaneous/ScrippsGroundingLine",
    processor=pooch.Unzip(),
)
groundingline: str = [file for file in shapefiles if file.endswith(".shp")][0]

# %%
# MEaSUREs Phase Map of Antarctic Ice Velocity
vel_file = f"{datafold}/Glaciology/MEaSUREs_PhaseBased_Velocity/antarctic_ice_vel_phase_map_v01"
vel: str = f"{vel_file}-vmag.nc"
try:
    assert os.path.exists(vel)
except AssertionError:
    # Note, download require a .netrc file containing 'machine urs.earthdata.nasa.gov login <uid> password <password>'
    # see https://nsidc.org/support/how/how-do-i-programmatically-access-data-spatial-temporal
    vel_file: str = pooch.retrieve(
        url="https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0754.001/1996.01.01/antarctic_ice_vel_phase_map_v01.nc",
        known_hash="fa0957618b8bd98099f4a419d7dc0e3a2c562d89e9791b4d0ed55e6017f52416",
        fname="antarctic_ice_vel_phase_map_v01.nc",
        path=f"{datafold}/Glaciology/MEaSUREs_PhaseBased_Velocity",
    )
    with pygmt.clib.Session() as lib:
        #! gmt grdmath ${vel_file}.nc?VX 2 POW ${vel_file}.nc?VY 2 POW POW 0.5 = ${vel_file}-vmag.nc
        lib.call_module(
            module="grdmath",
            args=f"{vel_file}.nc?VX 2 POW {vel_file}.nc?VY 2 POW POW 0.5 = {vel}",
        )

# %%
