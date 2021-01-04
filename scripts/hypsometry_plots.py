#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:57:37 2021

Hypsometry plots


@author: apj
"""


import pandas as pd
import geopandas as gpd
from rasterstats import zonal_stats
import rasterio as rio
import rasterio.mask
import matplotlib.pyplot as plt
import numpy as np
import glob
from shapely.validation import explain_validity
import os

# filepaths
fp_tiff = r'/Users/apj/Documents/_HY/Greenland/dem_diff/filled_dem/*.tif'
fp_outlines = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/*final.shp'
fp_diss_outline = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_final_edit.shp'
fp_exclude = r'/Users/apj/Documents/_HY/Greenland/dem_diff/filled_dem/glaciers_to_exclude_edit.csv'
fp_dem = r'/Users/apj/Documents/_HY/Greenland/DEM_masked/2016_dem_studyarea_3681x3295.tif'
fp_c = r'/Users/apj/Documents/_HY/Greenland/contour/2016_filled_contour.shp'


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


def areaDiff(outline, elevation_bin):
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


# function from pybob
def bin_data(bins, data2bin, bindata, mode='mean', nbinned=False):
    """
    Place data into bins based on a secondary dataset, and calculate statistics on them.
    :param bins: array-like structure indicating the bins into which data should be placed.
    :param data2bin: data that should be binned.
    :param bindata: secondary dataset that decides how data2bin should be binned. Should have same size/shape
        as data2bin.
    :param mode: How to calculate statistics of binned data. One of 'mean', 'median', 'std', 'max', or 'min'.
    :param nbinned: Return a second array, nbinned, with number of data points that fit into each bin.
        Default is False.
    :type bins: array-like
    :type data2bin: array-like
    :type bindata: array-like
    :type mode: str
    :type nbinned: bool
    :returns binned, nbinned: calculated, binned data with same size as bins input. If nbinned is True, returns a second
        array with the number of inputs for each bin.
    """
    assert mode in ['mean', 'median', 'std', 'max', 'min'], "mode not recognized: {}".format(mode)
    digitized = np.digitize(bindata, bins)
    binned = np.zeros(len(bins)) * np.nan
    if nbinned:  
        numbinned = np.zeros(len(bins))

    if mode == 'mean':
        for i, _ in enumerate(bins):
            binned[i] = np.nanmean(data2bin[np.logical_and(np.isfinite(bindata), digitized == i+1)])
            if nbinned:
                numbinned[i] = np.count_nonzero(np.logical_and(np.isfinite(data2bin), digitized == i+1))
    elif mode == 'median':
        for i, _ in enumerate(bins):
            binned[i] = np.nanmedian(data2bin[np.logical_and(np.isfinite(bindata), digitized == i+1)])
            if nbinned:
                numbinned[i] = np.count_nonzero(np.logical_and(np.isfinite(data2bin), digitized == i+1))
    elif mode == 'std':
        for i, _ in enumerate(bins):
            binned[i] = np.nanstd(data2bin[np.logical_and(np.isfinite(bindata), digitized == i+1)])
            if nbinned:
                numbinned[i] = np.count_nonzero(np.logical_and(np.isfinite(data2bin), digitized == i+1))
    elif mode == 'max':
        for i, _ in enumerate(bins):
            binned[i] = np.nanmax(data2bin[np.logical_and(np.isfinite(bindata), digitized == i+1)])
            if nbinned:
                numbinned[i] = np.count_nonzero(np.logical_and(np.isfinite(data2bin), digitized == i+1))
    elif mode == 'min':
        for i, _ in enumerate(bins):
            binned[i] = np.nanmin(data2bin[np.logical_and(np.isfinite(bindata), digitized == i+1)])
            if nbinned:
                numbinned[i] = np.count_nonzero(np.logical_and(np.isfinite(data2bin), digitized == i+1))
    else:
        raise ValueError('mode must be mean, median, std, max, or min')
    
    if nbinned:
        return np.array(binned), np.array(numbinned)
    else:
        return np.array(binned)

# function to extract data by mask
def glacierMask(fp_raster, features):
    with rasterio.open(fp_raster) as src:
        glac_mask, glac_out_transform = rasterio.mask.mask(src, shapes, crop=False)
        glac_nodata = src.nodata

    masked = glac_mask[0]
    masked[(masked == glac_nodata)] = np.nan    
    return masked

# function to get elevation bins
def getBins(array, bwidth):
    # get elev min and max
    elev_min = np.nanmin(array)
    elev_max = np.nanmax(array)
    # define elevation range
    erange = elev_max - elev_min
    
    min_el = elev_min - (elev_min % bwidth)
    max_el = elev_max + (bwidth - (elev_max % bwidth))
    bins = np.arange(min_el, max_el+1, bwidth)
    return bins

def excludeByID(excluded_list, in_gdf, id_column):
    # select all but ones in the list
    selection = in_gdf[~in_gdf[id_column].isin(excluded_list)]
    return selection

# read files
gdf = gpd.read_file(fp_diss_outline) # dissolved outlines
checkGeom(gdf)
contour = gpd.read_file(fp_c) # filled contours
exclude = pd.read_csv(fp_exclude) # glaciers to exclude
# rename columns
exclude = exclude.rename(columns={'Unnamed: 0': 'index', '0': 'ID'})
#remoce duplicates and drop extra columns
exclude = exclude.drop_duplicates(subset='ID')
exclude = exclude.drop(columns=['index'])
# to list
excludelist = exclude['ID'].tolist()

# glaciers to select
sel = ['RGI60-05.02328_1']
#sel = ['RGI60-05.02281', 'RGI60-05.02280_1', 'RGI60-05.02309', 'RGI60-05.02328_1', 'RGI60-05.01987_1', 'RGI60-05.02303', 'RGI60-05.02126'] # observed surges
#sel = ['RGI60-05.02281', 'RGI60-05.02280_1', 'RGI60-05.02328_1', 'RGI60-05.02309', 'RGI60-05.02108', 'RGI60-05.02303', 'RGI60-05.01920', 'RGI60-05.01987_1', 'RGI60-05.02213', 'RGI60-05.02126'] # observed and probable surge

# exclude and select certain glaciers 
gdf = excludeByID(excludelist, gdf, 'RGIId')
#gdf.to_file(fp_diss_outline_exc, driver='ESRI Shapefile')
gdf = gdf[gdf.RGIId.isin(sel)]
#gdf.to_file(fp_surging_outline, driver='ESRI Shapefile')

# read outline dataframes and assign them to dictionary where basename is the key to each dataframe
outlinedict = {} # empty dictionary
for f in glob.glob(fp_outlines):
    # get basename
    bname = os.path.basename(f)
    ol = gpd.read_file(f)
    ol = excludeByID(excludelist, ol, 'RGIId') # exclude certain glaciers
    ol = ol[ol.RGIId.isin(sel)] # subset dataframe
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

# exclude NoData from contour ranges
contsel = contour[contour['range_cont'] != '<NoData>']

# dissolve by contour range
contdis = contsel.dissolve(by='range_cont')
contdis = contdis.reset_index()

# read dem
with rio.open(fp_dem) as src:
    dem = src.read(1)
    dem_nodata = src.nodata

dem[dem == dem_nodata] = np.nan

# get geometries from selected polygons
shapes = gdf.geometry
    
# bins
bins = getBins(dem, 100)

# dataframe to store results
result = pd.DataFrame(bins, columns=['bins'])

# read tiff files to list
tifflist = []
for t in glob.glob(fp_tiff):
    tifflist.append(t)

# get elevation differences for each elevation bin
for tif in tifflist:    
    bname = os.path.basename(tif)
    # read ddem and mask with selected glaciers
    ddem = glacierMask(tif, shapes)
    
    # classify dem to bins
    digitized = np.digitize(dem, bins)
    
    # calculate average elevation difference per bin
    bin_means = bin_data(bins, ddem, dem, mode='mean', nbinned=False)        
    
    # parse column name
    colname = 'mu_dh_' + bname[0:12]
    
    # update results
    for i, _ in enumerate(bins):
        result.loc[result['bins'] == bins[i], colname] = bin_means[i]    
    
# update bins column to integer
result['bins'] = result['bins'].astype(int)

# add area change to new columns
# loop through dictionary keys and values
for x, y in outlinedict.items():
    # store the first four characters (the year) from the filename to variable
    year = x[:4]
    # dataframe to store results
    #result = pd.DataFrame(data=elev_bins, columns=['elev_bin'])
    # add column for results
    result[str(x[:4])+'Akm2'] = ""
    # loop through elevation bins and calculate area altitude difference for each bin
    for i in bins:
        i = i.astype(int)
        # selection by contour range before applying functions
        elev_bin = contdis[contdis['low_cont'] == i.astype(str)]
        # use function
        out = areaDiff(y, elev_bin)
        
        if out is None:
            out = 0
        # store result to dataframe
        result.loc[result['bins'] == i, str(x[:4])+'Akm2'] = out

# calculate area differences (e.g.2016 - 1953 so positive values show area increase and negative decrease)
result['dA53t85'] = result['1985Akm2'] - result['1953Akm2']
result['dA53t16'] = result['2016Akm2'] - result['1953Akm2']
result['dA85t16'] = result['2016Akm2'] - result['1985Akm2']

result = result.dropna(axis=0, how='any')

# figure output
fig_out = r'/Users/apj/Documents/_HY/Greenland/contour/figures/hypsometry_RGI60-0502328_surge_glaciers_filled_local_mean.png'

# create hypsometry and area altitude plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,10))

# hypsometry plot
line5385 = axes[0].plot(result['mu_dh_1953_to_1985'], result['bins'], marker='p', color='k', linewidth= 0.9, label='dh 1953 to 1985')
line532016 = axes[0].plot(result['mu_dh_1953_to_2016'], result['bins'], marker='v', color='b', linewidth=0.9, label='dh 1953 to 2016')
line852016 = axes[0].plot(result['mu_dh_1985_to_2016'], result['bins'], marker='s', color='g', linewidth=0.9, label='dh 1985 to 2016')
axes[0].set_ylim(min(result['bins'])-100, max(result['bins'])+100)
axes[0].set_yticks(np.arange(min(result['bins'])-100, max(result['bins'])+100, 100))
axes[0].axvline(0, color='grey', ls='--')
axes[0].set_ylabel('Elevation bin (m)')
axes[0].set_xlabel('Average elevation difference (m)')
axes[0].legend(loc=2)
axes[0].grid()

# area-altitude plot
area1953 = axes[1].plot(result['1953Akm2'], result['bins'], marker='s', color='k', linewidth= 0.9, label='1953')
area1985 = axes[1].plot(result['1985Akm2'], result['bins'], marker='^', color='#994C00', linewidth=0.9, label='1985')
area2016 = axes[1].plot(result['2016Akm2'], result['bins'], marker='o', color='#006633', linewidth=0.9, label='2016')
axes[1].set_ylim(min(result['bins'])-100, max(result['bins'])+100)
axes[1].set_yticks(np.arange(min(result['bins'])-100, max(result['bins'])+100, 100))
axes[1].axvline(0, color='grey', ls='--')
axes[1].set_ylabel('Elevation bin (m)')
axes[1].set_xlabel('Area altitude distribution ($km^2$)')
axes[1].legend(loc=1)
axes[1].grid()

# fig title
fig.suptitle('Glacier RGI60-05.02328', fontsize=20)
plt.tight_layout(pad=1.5)
plt.savefig(fig_out, dpi=600, format='png')

















