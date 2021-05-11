#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 13:21:13 2021


Hypsometry plot

To Do:
    select only glaciers where 'RGIId' matches with filtered dissolved outline 'RGIId'

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
fp_tiff = r'/Users/apj/Documents/_HY/Greenland/dem_diff/may/filled_ddem/*local*.tif'
fp_outlines = r'/Users/apj/Documents/_HY/Greenland/outlines/final_outlines_may/Dissolved_outlines.shp'
fp_diss_outline = r'/Users/apj/Documents/_HY/Greenland/outlines/final_outlines_may/Dissolved_outlines.shp'
fp_exclude = r'/Users/apj/Documents/_HY/Greenland/dem_diff/filled_ddem/glaciers_to_exclude_edit.csv'
fp_dem = r'/Users/apj/Documents/_HY/Greenland/kfs_dem/2016dem_wgs84_ellipsoid_utm_25m_studyarea_3681_3295.tif'

def checkGeom(geodataframe):
    """dit_
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
gdf = gpd.read_file(fp_diss_outline) # glacier outlines, the ones passing the filter in void fill
checkGeom(gdf)


# glaciers to select
#sel = ['RGI60-05.02328_1']
sel = ['RGI60-05.01920', 'RGI60-05.02108', 'RGI60-05.02213', 'RGI60-05.02328', 'RGI60-05.02280', 'RGI60-05.01987', 'RGI60-05.02216', 'RGI60-05.02297'] # active and probable glaciers
#sel = ['RGI60-05.02328', 'RGI60-05.02280', 'RGI60-05.01987', 'RGI60-05.02216', 'RGI60-05.02297'] # active glaciers
#sel = ['RGI60-05.01987_1', 'RGI60-05.01920', 'RGI60-05.02213'] # probable surge

# subset
gdf = gdf[gdf.RGIId.isin(sel)] 
# exclude and select certain glaciers 
#gdf = excludeByID(excludelist, gdf, 'RGIId')
# get id's to list, so same selection can be made for edited outlines
gdf_idlist = list(gdf['RGIId'])
#gdf.to_file(fp_diss_outline_exc, driver='ESRI Shapefile')


# read dem
with rio.open(fp_dem) as src:
    dem = src.read(1)
    dem_nodata = src.nodata

dem[dem == dem_nodata] = np.nan
  
# bins
bins = getBins(dem, 50)

# dataframe to store results
result = pd.DataFrame(bins, columns=['bins'])

# read tiff files to list
tifflist = []
for t in glob.glob(fp_tiff):
    tifflist.append(t)

# figure output
fig_out = r'/Users/apj/Documents/_HY/Greenland/dem_diff/may/figures/elev_diff_ind_glac_local.png'

# set figure outside loop so multiple plots are made to same figure
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,7))

# empty list for ids        
rgi_list = []

# counter
counter = 1

# loop through glaciers and plot elevation change in bins
for rgi in gdf.RGIId:
    
    # add rgi id to list
    rgi_list.append(rgi)
    
    # selection
    sel = gdf[gdf.RGIId == rgi]
    
    # get geometries from selected polygons
    shapes = sel.geometry
    
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
 
    # define line marker and color. simple way: define as many as selected glaciers to list and select by counter
    markerlist = ['1', '2', '.', 'o', 'v', 's', 'd', '+']
    colorlist = ['b', '#6b4c00', 'k', '#47ecff', '#EEE300', '#b96500', '#086800', '#ee33cc']
    
    #fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,6)) # this line is for testing plotting. Disable for actual run

    # hypsometry plot
    line5385 = axes[0].plot(result['bins'], result['mu_dh_1953_to_1985'], marker=markerlist[counter-1], color=colorlist[counter-1], linewidth= 0.9, markersize=3, label=rgi)
    line852016 = axes[1].plot(result['bins'], result['mu_dh_1985_to_2016'], marker=markerlist[counter-1], color=colorlist[counter-1], linewidth=0.9, markersize=3)
#    line532016 = axes[2].plot(result['bins'], result['mu_dh_1953_to_2016'], marker=markerlist[counter-1], color=colorlist[counter-1], linewidth=0.9, markersize=3)

    # fixed x axis
    axes[0].set_xlim(100, 2200)
    axes[0].set_xticks(np.arange(100,2200, 200) )
    axes[1].set_xlim(100, 2200)
    axes[1].set_xticks(np.arange(100,2200, 200) )
#    axes[2].set_xlim(100, 2200)
#    axes[2].set_xticks(np.arange(100,2200, 200) )
    
    # horizontal line at 0
    axes[0].axhline(0, color='grey', ls='--')
    axes[1].axhline(0, color='grey', ls='--')
#    axes[2].axhline(0, color='grey', ls='--')

    # titles
    axes[0].title.set_text('1953 - 1985')
    axes[1].title.set_text('1985 - 2016')
#    axes[2].title.set_text('1953 - 2016')
    # axis labels
#   axes[0].set_xlabel('Elevation bin', fontsize=15)
#    axes[1].set_xlabel('Elevation (m)', fontsize=15)
#    axes[2].set_xlabel('Elevation (m)', fontsize=15)
    # ylabel
#    axes[1].set_ylabel('Average elevation difference (m)', fontsize=15)
    # grids
#    axes[0].grid()
#    axes[1].grid()
#    axes[2].grid()
    # legend
    axes[0].legend(loc='upper right', bbox_to_anchor=(1.45, 1.01), fontsize = 13)
#    legend(loc='upper right', fontsize=13)

    # add 1 to counter
    counter += 1    

# add axes grid
axes[0].grid()
axes[1].grid()

# fig title
fig.text(0, 0.25, 'Average elevation difference (m)', rotation='vertical', fontsize=15)
fig.text(0.3, 0, 'Elevation (m)', fontsize=15)
fig.suptitle('Elevation changes on active and surge glaciers', fontsize=20, x=0.41)
plt.tight_layout(pad=3.5)
plt.savefig(fig_out, dpi=600, format='png')













