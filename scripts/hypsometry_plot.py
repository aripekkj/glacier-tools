#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 17:00:35 2020

Creating hypsometry plots from contour points with interpolated elevation values


@author: apj
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# filepath
fp = r'/Users/apj/Documents/_HY/Greenland/contour/53_contour_100_clip_pts_85and2016elev.shp'

# read file
gdf = gpd.read_file(fp)

# rename column
gdf = gdf.rename(columns={'RASTERVALU': '2016elev'})

# calculate elevation difference between Contour and interpolated elevation
gdf['85_ele_diff'] = gdf['85elev'] - gdf['CONTOUR']
gdf['2016_ele_diff'] = gdf['2016elev'] - gdf['CONTOUR']

# group by elevation contour and average values
grouped = gdf.groupby(['CONTOUR']).mean()

# create hypsometry plot
fig, ax = plt.subplots(figsize=(5,10))

line85 = ax.plot(grouped['85_ele_diff'], grouped.index, marker='p', color='k', linewidth= 0.9, label='1985')
line2016 = ax.plot(grouped['2016_ele_diff'], grouped.index, marker='v', color='b', linewidth=0.9, label='2016')
ax.set_ylim(min(grouped.index) - 100, max(grouped.index) + 100)
ax.yaxis.set_ticks(np.arange(min(grouped.index)-100, max(grouped.index)+100, 200))
ax.axvline(0, color='grey', ls='--')
ax.set_ylabel('Elevation (m)')
ax.set_xlabel('Average elevation difference to 1953 (m)')
ax.legend()
ax.grid()














