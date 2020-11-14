#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 17:00:35 2020

Creating hypsometry plots from contour points with interpolated elevation values


@author: apj
"""

import pandas as pd
import geopandas as gpd
from rasterstats import zonal_stats
import rasterio
import matplotlib.pyplot as plt
import numpy as np
import glob
from shapely.validation import explain_validity
import os

# filepaths
fp = r'/Users/apj/Documents/_HY/Greenland/contour/53_contour_100_clip_pts_85and2016elev.shp'
fp_c = r'/Users/apj/Documents/_HY/Greenland/contour/filled_contour.shp'
fp_outlines = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/*.shp'
fp_tiff = r'/Users/apj/Documents/_HY/Greenland/DEM_masked/coregistered/diff/*.tif'
fp_diss_outline = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm.shp'

def checkGeom(geodataframe):
    """
    Function to check validity of geometry. Returns message from shapely explain_validity if geometry is not 'Valid Geometry'

    Parameters
    ----------
    geodataframe : TYPE
        DESCRIPTION.

    Returns
    -------
    Message.

    """
    for geometry in geodataframe.geometry:
        if explain_validity(geometry) != 'Valid Geometry':
            print(explain_validity(geometry))

# read file
gdf = gpd.read_file(fp_diss_outline)
checkGeom(gdf)
contour = gpd.read_file(fp_c)

# read tiff files to list
tifflist = []
for t in glob.glob(fp_tiff):
    tifflist.append(t)

# read outline dataframes and assign them to dictionary where basename is the key to each dataframe
outlinedict = {} # empty dictionary
for f in glob.glob(fp_outlines):
    # get basename
    bname = os.path.basename(f)
    ol = gpd.read_file(f)
    # check geometry validity before creating dictionary
    for geom in ol.geometry:
        if explain_validity(geom) != 'Valid Geometry':
            print(bname + ' Geometry has invalid parts')
            print(explain_validity(geom))
    # add dataframe to dictionary with basename as key
    outlinedict[bname] = ol

# check contour columns and elevation ranges
contour.columns
contour['range_cont'].unique()

# select only positive contour ranges with interval of 100
contsel = contour[contour['range_cont'] != '-100 - 0']
contsel = contsel[contsel['range_cont'] != '0 - 400']

# dissolve by contour range
contdis = contsel.dissolve(by='range_cont')
contdis = contdis.reset_index()

# set contour range, for example 0 - 100
#contour_range = '0 - 100'

# select contour range before applying functions
#elev_bin = contdis[contdis['range_cont'] == contour_range]


def areaDiff(outline, elevation_bin, contour_range):
    """
    Function to calculate area in an elevation bin

    Parameters
    ----------
    outline : Polygon
        Polygon containing outlines.
    elevation_bin : Polygon
        Polygon containing elevation ranges
    contour_range : String
        Elevation range to be selected

    Returns
    -------
    elev_range_area_sum : float
        Sum of areas from outline polygon inside the elevation bin

    """
    # clip outlines by selected elevation range
    outline_elev_range = gpd.clip(outline, elevation_bin, keep_geom_type=(True))
    # check that clipped dataframe is not empty
    if outline_elev_range.empty == True:
        return
    # compute area in km2
    elev_range_area = outline_elev_range.geometry.area / 1000000
    # sum areas
    elev_range_area_sum = elev_range_area.sum()
    return elev_range_area_sum


def areaAltDiff(outline, elevation_bin, contour_range, diff_raster):
    """
    Function to calculate average elevation change in elevation bin

    Parameters
    ----------
    outline : Polygon
        Polygon containing the outlines of the area.
    elevation_bin : Polygon
        Polygon of the elevation range.
    diff_raster : Raster
        Difference raster of two DEM surfaces.

    Returns
    -------
    
    zonal_mean : GeoDataFrame
        Zonal mean of input raster in the elevation bin

    """
    # clip outlines by selected elevation range
    outline_elev_range = gpd.clip(outline, elevation_bin, keep_geom_type=(True))
    # check that clipped dataframe is not empty
    if outline_elev_range.empty == True:
        print('No features in DataFrame after clipping')
        return
    
    # calculate average elevation change in the elevation bin
    with rasterio.open(diff_raster) as src:
        array = src.read(1)
        affine = src.transform
    # compute zonal mean and get output as geojson
    zon_mean = zonal_stats(outline_elev_range, array, affine=affine, stats='mean', all_touched=True, nodata=src.nodata, geojson_out=True)
    # turn geojson to GeoDataFrame
    zonal_mean = gpd.GeoDataFrame.from_features(zon_mean)
    # assign elevation bin to column
    zonal_mean['elev_bin'] = contour_range
    # group by elev_bin
    grouped = zonal_mean.groupby(['elev_bin']).mean()
    # reset index
    grouped = grouped.reset_index()
    # keep only columns 'mean' , 'elev_bin'
    grouped = grouped[['elev_bin', 'mean']]

    return grouped

# create list of elevation bins
elev_bins = list(contdis['range_cont'].unique())

# empty list for result dataframes
#resultlist = []

resultdf = pd.DataFrame(data=elev_bins, columns=['elev_bin'])
###### Get total area for each elevation bin ######
# loop through dictionary keys and values
for x, y in outlinedict.items():
    # store the first four characters (the year) from the filename to variable
    year = x[:4]
    # dataframe to store results
    result = pd.DataFrame(data=elev_bins, columns=['elev_bin'])
    # add column for results
    result[str(x[:4])+'Akm2'] = ""
    # loop through elevation bins and calculate area altitude difference for each bin
    for i in elev_bins:
        # selection by contour range before applying functions
        elev_bin = contdis[contdis['range_cont'] == i]
        # use function
        out = areaDiff(y, elev_bin, i)
        # store result to dataframe
        result.loc[result['elev_bin'] == i, str(x[:4])+'Akm2'] = out
        
    if len(resultdf) == 0:
        resultdf = result
    else:
        resultdf = resultdf.merge(result, how='inner', on='elev_bin')
    
    # add result dataframe to list
    #resultlist.append(result)

zonmean_result = pd.DataFrame(data=elev_bins, columns=['elev_bin'])
##### Get zonal mean difference ######
for tiffile in tifflist:
    tifname = os.path.basename(tiffile[:-4])
    zonmean = pd.DataFrame(data=elev_bins, columns=['elev_bin'])
    # add column for results
    zonmean['mu_dh_'+tifname] = ""
    for i in elev_bins:
        # selection by contour range before applying functions
        elev_bin = contdis[contdis['range_cont'] == i]
        output = areaAltDiff(gdf, elev_bin, i, tiffile)
        # store output to df
        if output is None:
            zonmean.loc[zonmean['elev_bin'] == i, 'mu_dh_'+tifname] = output
        else:    
            zonmean.loc[zonmean['elev_bin'] == i, 'mu_dh_'+tifname] = output['mean'][0]
    
    if len(zonmean_result) == 0:
        zonmean_result = zonmean
    else:
        zonmean_result = zonmean_result.merge(zonmean, how='inner', on='elev_bin')
    print('End of loop')
    # update column name
    #zonmean_result = zonmean_result.rename(columns={'mean': 'mu_dh_' + tifname})    


# merge area and zonal mean results
merged = pd.merge(resultdf, zonmean_result, how='inner', on='elev_bin')
merged.columns
# create column for high contour
merged['h_cont'] = merged['elev_bin'].str.split(expand=True)[2].astype(str).astype(int)

# sort ascending
merged = merged.sort_values(by='h_cont', ascending=True)
merged

# calculate area differences (e.g.2016 - 1953 so positive values show area increase and negative decrease)
merged['dA53t85'] = merged['1985Akm2'] - merged['1953Akm2']
merged['dA53t16'] = merged['2016Akm2'] - merged['1953Akm2']
merged['dA85t16'] = merged['2016Akm2'] - merged['1985Akm2']


# create hypsometry and area altitude plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,10))

# hypsometry plot
line85 = axes[0].plot(merged['mu_dh_53to85'], merged['h_cont'], marker='p', color='k', linewidth= 0.9, label='dh 1953 to 1985')
line2016 = axes[0].plot(merged['mu_dh_53to2016'], merged['h_cont'], marker='v', color='b', linewidth=0.9, label='dh 1953 to 2016')
axes[0].set_ylim(min(merged['h_cont'])+100, max(merged['h_cont'])+200)
axes[0].set_yticks(np.arange(min(merged['h_cont'])+100, max(merged['h_cont'])+200, 100))
axes[0].axvline(0, color='grey', ls='--')
axes[0].set_ylabel('Elevation (m)')
axes[0].set_xlabel('Average elevation difference to 1953 (m)')
axes[0].legend()
axes[0].grid()

# area-altitude plot
line5385 = axes[1].plot(merged['dA53t85'], merged['h_cont'], marker='s', color='k', linewidth= 0.9, label='area change 1953 to 1985')
line532016 = axes[1].plot(merged['dA53t16'], merged['h_cont'], marker='^', color='r', linewidth=0.9, label='area change 1953 to 2016')
axes[1].set_ylim(min(merged['h_cont'])+100, max(merged['h_cont'])+200)
axes[1].set_yticks(np.arange(min(merged['h_cont'])+100, max(merged['h_cont'])+200, 100))
axes[1].axvline(0, color='grey', ls='--')
axes[1].set_ylabel('Elevation (m)')
axes[1].set_xlabel('Area difference to 1953 (km^2)')
axes[1].legend()
axes[1].grid()











