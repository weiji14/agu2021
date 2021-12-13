# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.13.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # **3D point cloud of Siple Coast active subglacial lakes**
#
# Making a LiDAR point cloud of Antarctica for an 2021 AGU poster.
# Will be uploaded as a 3D model to [Sketchfab](https://sketchfab.com).

# %%
import cmcrameri.cm as cmc
import pandas as pd
import pygmt
import pyvista as pv


# %% [markdown]
# ## Load ICESat-2 point cloud over Siple Coast

# %%
# Load Siple Coast pre-processed ATL11 point cloud data
# This 4.6GB file was processed using the atlxi_dhdt.ipynb script from
# https://github.com/weiji14/deepicedrain/pull/329/commits/9073678a2adb42fcc863afec22957c879988fa3d
df_dhdt: pd.DataFrame = pd.read_parquet(path="df_dhdt_siple_coast.parquet")

# %%
# Filter point cloud to those with rate of elevation change (dhdt)
# that is -0.12 m/yr < dhdt_slope > 0.12 m/yr
_df_dhdt = df_dhdt[abs(df_dhdt.dhdt_slope) > 0.12]  # 0.105
points: pd.DataFrame = (
    _df_dhdt[["x", "y", "h_corr_11", "dhdt_slope"]]
    .rename(columns=dict(h_corr_11="z"))
    .dropna()
)

# Add some vertical exaggeration (x100) to z-axis
points["z"] *= 100
points

# %% [markdown]
# ## Wrap ICESat-2 data in PyVista format

# %%
# Create XYZ point cloud in pyvista format
cloud = pv.PolyData(var_inp=points[["x", "y", "z"]].to_numpy())

# %%
# Add dhdt_slope as an attribute to the XYZ point cloud
cloud.point_data["dhdt_slope"] = points.dhdt_slope
cloud

# %% [markdown]
# ## Render the point cloud!

# %%
# Make a 3D PyVista plot of the XYZ point cloud, colored by dhdt_slope
p = pv.Plotter()
p.add_points(points=cloud, clim=[-2.5, 2.5], cmap=cmc.vik_r)
p.export_gltf(filename="siple_coast_point_cloud.gltf")
p.show(cpos="xy")


# %%
# https://superuser.com/questions/281573/what-are-the-best-options-to-use-when-compressing-files-using-7-zip/742034#742034
# !7z a -t7z -m0=lzma2 -mx=9 -mfb=64 -md=32m -ms=on siple_coast_point_cloud.gltf.7z siple_coast_point_cloud.gltf

# %%
# Quick plot
# cloud.plot(cpos="xy", render_points_as_spheres=True, clim=[-2.5, 2.5], cmap=cmc.vik_r)
