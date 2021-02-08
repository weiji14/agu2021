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
import re

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
    try:
        # Run processing of velocity magnitude grid locally
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
    except:
        # Just download pre-processed velocity magnitude grid from GitHub
        vel: str = pooch.retrieve(
            url="https://github.com/weiji14/nzasc2021/releases/download/v0.0.0/antarctic_ice_vel_phase_map_v01-vmag.nc",
            known_hash="ed6393275d8d8475c2162a838d6b9220cd529d28a2b5d674a6bf6dbda4971049",
            fname="antarctic_ice_vel_phase_map_v01-vmag.nc",
            path=f"{datafold}/Glaciology/MEaSUREs_PhaseBased_Velocity",
        )

# %%
# DeepIceDrain active subglacial lake outlines
lakes = "https://raw.githubusercontent.com/weiji14/deepicedrain/v0.4.0/antarctic_subglacial_lakes_3031.geojson"
lakes = "https://github.com/weiji14/deepicedrain/blob/6bbd5831b2fa9bfaf69732aaad6fa6822e02a8d0/antarctic_subglacial_lakes_3031.geojson"


# %% [markdown]
# ### Make color maps
#
# For MOA and ice velocity and vertical elevation trend (dhdt).
# Scientific color maps are from http://doi.org/10.5281/zenodo.1243862

# %%
pygmt.makecpt(
    series=[15000, 17000, 1],
    cmap="grayC",
    continuous=True,
    output="cmap_moa.cpt",
    reverse=True,
)
with pygmt.config(COLOR_FOREGROUND="240/249/33", COLOR_BACKGROUND="13/8/135"):
    pygmt.makecpt(series=[0, 800, 1], cmap="batlow", output="cmap_vel.cpt")
pygmt.makecpt(
    cmap="berlin",
    series=[-3, 3, 1],
    reverse=True,
    continuous=True,
    output="cmap_dhdt.cpt",
)

# %%

# %% [markdown]
# # Figure of Siple Coast active subglacial lakes

# %%
# We're making this a specific height
figheight = 115  # in mm

# Region in PS71 for main part of the figure
sip_xl, sip_xh, sip_yl, sip_yh = sipreg = [-800_000, 25_000, -1_000_000, -400_000]

# Calculate the figure width and map scale
figwidth = figheight * (sip_xh - sip_xl) / (sip_yh - sip_yl)
sipratio = (sip_yh - sip_yl) / (figheight / 1000)

# Make a GMT region string and projection strings in both PS71 and Lon/Lat
sipreg = [sip_xl, sip_xh, sip_yl, sip_yh]
sipproj = f"x1:{sipratio}"
sipproj_ll = f"s0/-90/-71/1:{sipratio}"

# %%
# Initialize figure and plot MOA as the basemap with ticks every 200 km in xy directions
fig = pygmt.Figure()
with pygmt.config(MAP_FRAME_TYPE="inside"):
    fig.basemap(
        region=sipreg, projection=sipproj, frame=["nwse", "xf200000", "yf200000"]
    )
    fig.grdimage(grid=moa, cmap="cmap_moa.cpt", nan_transparent=True)

# Plot graticules overtop, every 2° latitude and 15° longitude
with pygmt.config(
    MAP_ANNOT_OFFSET_PRIMARY="-2p",
    MAP_FRAME_TYPE="inside",
    MAP_ANNOT_OBLIQUE=0,
    FONT_ANNOT_PRIMARY="8p,grey",
    MAP_GRID_PEN_PRIMARY="grey",
    MAP_TICK_LENGTH_PRIMARY="-10p",
    MAP_TICK_PEN_PRIMARY="thinnest,grey",
    FORMAT_GEO_MAP="dddF",
    MAP_POLAR_CAP="88/90",  # less longitude graticules at >88°S
):
    fig.basemap(
        projection=sipproj_ll, region=sipreg, frame=["NSWE", "xa15g15", "ya2g2"]
    )

# Plot the grounding line in white
fig.plot(data=groundingline, region=sipreg, projection=sipproj, pen="0.15p,white")


# %%
# Overlay ice velocity with 70% transparency
fig.grdimage(grid=vel, cmap="cmap_vel.cpt", transparency=70, nan_transparent=True)
# Overlay dhdt with 30% transparency
# pygmt.makecpt(cmap="berlin", series=[-1.0, 1.0, 0.25], continuous=True, reverse=True)
# fig.grdimage(
#     grid="ATLXI/ds_grid_dhdt_siple_coast.nc",
#     cmap=True,
#     # cmap="cmap_dhdt.cpt",
#     transparency=30,
#     nan_transparent=True,
# )
fig.show()


# %%
# Plot lakes in PS71 as blobs (red for draining, blue for filling)
# TODO refactor after pygmt/geopandas integration is done,
# see https://github.com/GenericMappingTools/pygmt/issues/608
with pygmt.helpers.GMTTempFile(suffix=".gmt") as tmpfile:
    os.remove(path=tmpfile.name)
    gpd.read_file(lakes).to_file(tmpfile.name, driver="OGR_GMT")
    fig.plot(
        data=tmpfile.name,  # "antarctic_subglacial_lakes_3031.geojson"
        pen="thinnest,yellow,-",
        cmap="cmap_dhdt.cpt",
        color="+z",
        close=True,
        # transparency=30,
        a="Z=inner_dhdt",
    )

# %%
# Siple Coast placename labels
gdf = gpd.read_file("antarctic_subglacial_lakes_3031.gmt")
with open("place_labels_siple_coast.tsv", mode="w") as file:
    # Ice Streams A to E
    font = "7p,Helvetica-Narrow-Oblique,white"
    print(f"-320000\t-440000\t-65\t{font}\tCM\tMercer Ice Stream", file=file)
    print(f"-385000\t-555000\t-5\t{font}\tCM\tWhillans Ice Stream", file=file)
    print(f"-470000\t-450000\t-55\t{font}\tCM\tVan der Veen Ice Stream", file=file)
    print(f"-550000\t-625000\t-40\t{font}\tCM\tKamb Ice Stream", file=file)
    print(f"-700000\t-700000\t-45\t{font}\tCM\tBindschadler Ice Stream", file=file)
    print(f"-700000\t-850000\t-37\t{font}\tCM\tMacAyeal Ice Stream", file=file)

    # Ice ridges, rises and domes
    print(f"-370000\t-480000\t-30\t{font}\tCM\tConway Ice Ridge", file=file)
    print(f"-400000\t-600000\t0\t{font}\tCM\tEngelhardt Ice Ridge", file=file)
    print(f"-100000\t-750000\t-45\t{font}\tCM\tCrary Ice Rise", file=file)
    print(f"-450000\t-780000\t-35\t{font}\tCM\tSiple Dome", file=file)
    print(f"-650000\t-800000\t-25\t{font}\tCM\tShabtaie Ice Ridge", file=file)
    print(f"-650000\t-950000\t-15\t{font}\tCM\tHarrison Ice Ridge", file=file)

    # Abbreviated lakes
    abbrev_dict: dict = {
        "Subglacial": "S",
        "Lake": "L",
        "Conway": "C",
        "Engelhardt": "E",
        "Kamb": "K",
        "MacAyeal": "Mac",
        "Mercer": "M",
        # "Recovery", "R",
        # "Slessor","S"
        "Whillans": "W",
    }
    pattern = re.compile(pattern=rf"\b({'|'.join(abbrev_dict.keys())})\b")

    for idx, row in gdf.dissolve(by="lake_name", as_index=False).iterrows():
        x, y = row.geometry.centroid.xy

        label = pattern.sub(
            repl=lambda name: abbrev_dict[name.group()], string=row.lake_name
        ).replace(" ", "")
        justify = (
            "BL"
            if label in ("L12", "L78", "SLW", "W7", "WX", "WXI")
            else "TC"
            if ("*" in label or label == "SLM")
            else "TR"
        )

        print(f"{x[0]:.0f}\t{y[0]:.0f}\t0\t6p,white\t{justify}\t{label}", file=file)

# %%
# Plot labels for Siple Coast ice streams, active subglacial lakes, etc
fig.text(
    textfiles="place_labels_siple_coast.tsv",
    angle=True,
    font=True,
    justify=True,
    offset="j0.12c",
    # frame=["WsNe", "af10000g50000"],
)
fig.show()

# %%
# Plot the color bar once with a transparent box, then again with no box and no transparency
with pygmt.config(
    FONT_ANNOT_PRIMARY="6p,white",
    FONT_LABEL="6p,white",
    MAP_ANNOT_OFFSET_PRIMARY="2p",
    MAP_TICK_PEN_PRIMARY="0.25p,white",
    MAP_TICK_LENGTH_PRIMARY="3p",
    MAP_FRAME_PEN="0.5p,white",
    MAP_LABEL_OFFSET="4p",
):
    colorbar_kwargs = dict(
        cmap="cmap_dhdt.cpt",
        position="jBR+jBR+w1.6c/0.18c+o1.2c/0.3c+v+e",
        frame=["xaf", 'y+l"dhdt (m/yr)"'],
    )
    fig.colorbar(
        box="+gblack+c-9p/3p",
        # box = '+gblack+p0.5p,black+c3p'
        transparency=70,
        **colorbar_kwargs,
    )
    fig.colorbar(**colorbar_kwargs)

    # Add a scalebar
    fig.basemap(
        projection=sipproj_ll, region=sipreg, map_scale="jBR+o2.2c/0.3c+w50k+uk+f"
    )

fig.show()

# %%
# Make insets of Antarctica
# Workaround until https://github.com/GenericMappingTools/pygmt/pull/788 is merged
antwidth = 3  # width of inset in cm
with pygmt.clib.Session() as lib:
    lib.call_module(module="inset", args=f"begin -DjTR+w{antwidth}c")

    # Plot the inset map
    fig.basemap(region=vel, projection=f"S0/-90/71/{antwidth}c", frame="+n")
    fig.coast(area_thresh="+a", land="white")  # ice shelf in white
    fig.coast(area_thresh="+ag", land="gray")  # grounded ice in gray
    fig.plot(
        projection=f"X{antwidth}c",
        x=[sip_xl, sip_xl, sip_xh, sip_xh, sip_xl],
        y=[sip_yl, sip_yh, sip_yh, sip_yl, sip_yl],
        pen="1p,black",  # map location in black
    )

    lib.call_module(module="inset", args="end")

fig.show()

# %%
# Save the figure
fig.savefig(fname="siple_coast_lakes.pdf")
fig.savefig(fname="siple_coast_lakes.png", dpi=1200)


# %%

# %% [markdown]
# # Figure of Antarctic active subglacial lake map

# %%
# We're making this a specific height
figheight = 115  # in mm

# Region in PS71 for main part of the figure
ais_xl, ais_xh, ais_yl, ais_yh = aisreg = [-2700000, 2800000, -2200000, 2300000]

# Calculate the figure width and map scale
figwidth = figheight * (ais_xh - ais_xl) / (ais_yh - ais_yl)
aisratio = (ais_yh - ais_yl) / (figheight / 1000)

# Make a GMT region string and projection strings in both PS71 and Lon/Lat
aisreg = [ais_xl, ais_xh, ais_yl, ais_yh]
aisproj = "x1:" + str(aisratio)
aisproj_ll = "s0/-90/-71/1:" + str(aisratio)


# %%
# Initialize figure and plot MOA as the base map with ticks every 200 km both directions
fig = pygmt.Figure()
with pygmt.config(MAP_FRAME_TYPE="inside"):
    # fig.coast(
    #     projection=aisproj_ll, region=aisreg, land="lightblue", water="royalblue2"
    # )
    fig.coast(region=aisreg, projection=aisproj_ll, resolution="c", water=True)
    fig.grdimage(
        grid="@earth_relief_03m",
        projection=aisproj_ll,
        region=aisreg,
        cmap="oleron",
        shading=True,
    )
    fig.coast(Q=True)  # end water clip path
    fig.basemap(
        projection=aisproj, region=aisreg, frame=["nwse", "xf200000", "yf200000"]
    )
    fig.grdimage(grid=moa, cmap="cmap_moa.cpt", nan_transparent=True)

# Plot graticules overtop, every 10° latitude and 45° longitude
with pygmt.config(
    MAP_ANNOT_OFFSET_PRIMARY="-2p",
    MAP_FRAME_TYPE="inside",
    MAP_ANNOT_OBLIQUE=0,
    FONT_ANNOT_PRIMARY="8p,grey",
    MAP_GRID_PEN_PRIMARY="grey",
    MAP_TICK_LENGTH_PRIMARY="-10p",
    MAP_TICK_PEN_PRIMARY="thinnest,grey",
    FORMAT_GEO_MAP="dddF",
    MAP_POLAR_CAP="88/90",  # less longitude graticules at >88°S
):
    fig.basemap(
        projection=aisproj_ll, region=aisreg, frame=["NSWE", "xa45g45", "ya10g10"]
    )

# Plot the grounding line in white
fig.plot(data=groundingline, region=aisreg, projection=aisproj, pen="0.15p,white")

# Plot bounding box of Siple Coast study area
fig.plot(
    x=[sip_xl, sip_xl, sip_xh, sip_xh, sip_xl],
    y=[sip_yl, sip_yh, sip_yh, sip_yl, sip_yl],
    pen="0.5p,black",  # map location in black
    transparency=50,
)

# %%
# Overlay ice velocity with 70% transparency
fig.grdimage(grid=vel, cmap="cmap_vel.cpt", transparency=70, nan_transparent=True)
fig.show()

# %%
# Antarctica placename labels
with open("place_labels_antarctica.tsv", mode="w") as file:
    font = "10p,Helvetica-Narrow,white"
    print(f"-1070000\t-360000\t0\t{font}\tWest", file=file)
    print(f"-1070000\t-460000\t0\t{font}\tAntarctica", file=file)
    print(f"950000\t600000\t0\t{font}\tEast", file=file)
    print(f"950000\t500000\t0\t{font}\tAntarctica", file=file)

    font = "6p,Helvetica-Narrow-Oblique,white"
    print(f"0\t-1000000\t0\t{font}\tRoss", file=file)
    print(f"0\t-1060000\t0\t{font}\tIce Shelf", file=file)
    print(f"-1000000\t700000\t0\t{font}\tRonne-Filchner", file=file)
    print(f"-1000000\t640000\t0\t{font}\tIce Shelf", file=file)
    print(f"2100000\t760000\t0\t{font}\tAmery", file=file)
    print(f"2100000\t700000\t0\t{font}\tIce Shelf", file=file)

    font = "4p,Helvetica-Narrow-Oblique,white"
    print(f"-400000\t1150000\t45\t{font}\tSlessor Glacier", file=file)
    print(f"-90000\t850000\t0\t{font}\tRecovery Glacier", file=file)
    print(f"-340000\t200000\t-30\t{font}\tFoundation Ice Stream", file=file)
    print(f"-810000\t95000\t-60\t{font}\tInstitute Ice Stream", file=file)
    # print(f"-1500000\t50000\t20\t{font}\tRutford Ice Stream", file=file)
    print(f"-1300000\t-300000\t30\t{font}\tThwaites Glacier", file=file)
    print(f"-500000\t-700000\t40\t{font}\tSiple Coast", file=file)
    print(f"500000\t-440000\t0\t{font}\tNimrod 2", file=file)
    print(f"550000\t-750000\t35\t{font}\tByrd Glacier", file=file)
    print(f"700000\t-1500000\t5\t{font}\tDavid Glacier", file=file)
    print(f"870000\t-1700000\t0\t{font}\tCook E2", file=file)
    # print(f"2200000\t-850000\t-70\t{font}\tTotten Glacier", file=file)
    print(f"1400000\t650000\t15\t{font}\tLambert Glacier", file=file)

# %%
# Plot labels for Antarctic ice shelves, ice streams, etc
fig.text(
    region=aisreg,
    projection=aisproj,
    textfiles="place_labels_antarctica.tsv",
    angle=True,
    font=True,
    # frame=["WsNe", "af100000g500000"],
)

# %%
# Plot lakes in PS71 as cyan blobs with 60% transparency
# TODO refactor after pygmt/geopandas integration is done,
# see https://github.com/GenericMappingTools/pygmt/issues/608
with pygmt.helpers.GMTTempFile(suffix=".gmt") as tmpfile:
    os.remove(path=tmpfile.name)
    gpd.read_file(lakes).to_file(tmpfile.name, driver="OGR_GMT")
    fig.plot(
        data=tmpfile.name,  # "antarctic_subglacial_lakes_3031.geojson"
        pen="0.5p,cyan",
        color="cyan",
        transparency=60,
    )
fig.show()

# %%
# Add the finishing touches by putting a panel label in the bottom left corner,
# and add a legend with a 70% transparent box behind it (so you can see it
# clearly) that includes:
# - a cyan ellipse the cyan lakes
# - a color bar for the ice surface velocity
barwidth = 0.3 * figwidth / 10  # legend width is 30% of the map width, in cm

# Position string for feeding to the colorbar call. A bit coded, but the pyGMT
# documents are great: https://www.pygmt.org/dev/api/index.html
legend_pos: str = f"jBL+jBL+w{barwidth}c+o0.2c/0.2c"

# Plot the color bar once with a transparent box, then again with no box and no transparency
with pygmt.config(
    FONT_ANNOT_PRIMARY="8p,white",
    FONT_LABEL="8p,white",
    MAP_ANNOT_OFFSET_PRIMARY="2p",
    MAP_TICK_PEN_PRIMARY="0.5p,white",
    MAP_TICK_LENGTH_PRIMARY="3p",
    MAP_FRAME_PEN="0.5p,white",
    MAP_LABEL_OFFSET="4p",
):
    fig.legend(
        spec="legend.txt", position=legend_pos, box="+gblack+c1.2p", transparency="70"
    )
    fig.legend(spec="legend.txt", position=legend_pos)

fig.show()

# %%
# Save the figure
fig.savefig(fname="antarctica_lakes.pdf")
fig.savefig(fname="antarctica_lakes.png", dpi=600)

# %%
