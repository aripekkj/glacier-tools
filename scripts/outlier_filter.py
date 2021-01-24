#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:35:03 2020

Outlier removal function from pybob package
    https://pybob.readthedocs.io/en/stable/modules/pybob.GeoImg.html

Script saves results of 3-sigma filter and nmad filter

@author: apj
"""

import fiona
import numpy as np
import rasterio as rio
import rasterio.mask
import matplotlib.pyplot as plt
import geopandas as gpd
import os

# filepaths
fp_dem = r'/Users/apj/Documents/_HY/Greenland/kfs_dem/2016dem_wgs84_ellipsoid_utm_25m_studyarea_edit.tif' # dem where to derive elevation bins
#fp_ddem = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_01122020/test/1985_2016_02037.tif'
fp_ddem = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/edited_studyarea/1953_to_1985_ddem.tif' # ddem to filter
out_dir = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filtered/'
fp_glac = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_final.shp'

# read shapefile geometries
with fiona.open(fp_glac, "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

# read raster files with rasterio
with rio.open(fp_dem) as src:
    dem = src.read(1)
    p = src.profile
    dem_nodata = src.nodata

#with rio.open(fp_ddem) as src:
#    ddem = src.read(1)
#    p_ddem = src.profile
#    ddem_nodata = src.nodata

with rasterio.open(fp_ddem) as src:
    #ddem = src.read(1)
    out_image, out_transform = rasterio.mask.mask(src, shapes, crop=False) # mask with polygons
    ddem_nodata = src.nodata
    p_ddem = src.profile
    out_meta = src.meta

with rasterio.open(fp_ddem) as src:
    inverse_out_image, inverse_out_transform = rasterio.mask.mask(src, shapes, crop=False, invert=True)
    inverse_out_meta = src.meta


# function to get elevation bins
def getBins(array, bwidth):
    # get elev min and max
    elev_min = np.nanmin(array)
    elev_max = np.nanmax(array)
    # define elevation range
    
    min_el = elev_min - (elev_min % bwidth)
    max_el = elev_max + (bwidth - (elev_max % bwidth))
    bins = np.arange(min_el, max_el+1, bwidth)
    return bins

#plt.imshow(out_image[0], cmap='BuPu')
#plt.show()

ddem = out_image[0]
ddem_outside_glaciers = inverse_out_image[0]

# change no data to nan
ddem[(ddem == ddem_nodata)] = np.nan
dem[(dem == dem_nodata)] = np.nan
ddem_outside_glaciers[(ddem_outside_glaciers == ddem_nodata)] = np.nan

# create array of bins
#bins = np.arange(0, np.nanmax(dem), 100) #np.max(dem)
bins = getBins(dem, 100)

""" testing
for i, _ in enumerate(bins):
    print(i)
    
i = 0
new_dem = np.zeros(dem.shape)
digitized = np.digitize(dem, bins) # elevation bins
this_bindata = dem[digitized == i]
new_dem[digitized == i] = this_bindata
plt.imshow(new_dem, cmap='BuPu')
plt.show()



# new dem
new_dem = np.zeros(dem.shape) # array with zeros of the size of dem array
digitized = np.digitize(dem, bins) # elevation bins


#Finding out what happens inside the function
nsig=3
for i, _ in enumerate(bins):
    print(i)
    if i == 1:
        break
    this_bindata = dem[digitized == i]
    nout = 1
    old_mean = np.nanmean(this_bindata)
    old_std = np.nanstd(this_bindata)
    while nout >= 1:
        thresh_up = old_mean + nsig * old_std
        thresh_dn = old_mean - nsig * old_std
        print('bin ' + str(i))
        print('Upper threshold: ' + str(thresh_up))
        print('Lower threshold: ' + str(thresh_dn))


        isout = np.logical_or(this_bindata > thresh_up, this_bindata < thresh_dn)
        nout = np.count_nonzero(isout)
        this_bindata[isout] = np.nan
        
        old_mean = np.nanmean(this_bindata)
        old_std = np.nanstd(this_bindata)
    new_dem[digitized == i] = this_bindata
"""


def outlier_removal(bins, DEM, dDEM, nsig=3):
    """
    Iteratively remove outliers in an elevation bin using a 3-sigma filter.
    :param bins: lower bound of elevation bins to use
    :param DEM: DEM to determine grouping for outlier values
    :param dDEM: elevation differences to filter outliers from
    :param nsig: number of standard deviations before a value is considered an outlier.
    :type bins: array-like
    :type DEM: array-like
    :type dDEM: array-like
    :type nsig: float
    :returns new_ddem: ddem with outliers removed (set to NaN)
    """
    new_ddem = np.zeros(dDEM.shape)
    digitized = np.digitize(DEM, bins)
    for i, _ in enumerate(bins):
        # Option to exclude elevation bins from outlier removal
        #if i == 0:
         #   print('Excluding outlier filtering on elevations lower than ' + str(bins[0]))
          #  this_bindata = dDEM[digitized == i]
           # new_ddem[digitized == i] = this_bindata
           # continue
        this_bindata = dDEM[digitized == i]
        nout = 1
        old_mean = np.nanmean(this_bindata)
        old_std = np.nanstd(this_bindata)
        while nout >= 1:
            thresh_up = old_mean + nsig * old_std
            thresh_dn = old_mean - nsig * old_std
            #print('bin ' + str(bins[i]))
            #print('Upper threshold: ' + str(thresh_up))
            #print('Lower threshold: ' + str(thresh_dn))

            isout = np.logical_or(this_bindata > thresh_up, this_bindata < thresh_dn)
            nout = np.count_nonzero(isout)
            this_bindata[isout] = np.nan

            old_mean = np.nanmean(this_bindata)
            old_std = np.nanstd(this_bindata)
        new_ddem[digitized == i] = this_bindata
        plt.imshow(new_ddem, cmap='BuPu')
        plt.show()

    return new_ddem

def nmad(data, nfact=1.4826):
    """
    Calculate the normalized median absolute deviation (NMAD) of an array.
    :param data: input data
    :param nfact: normalization factor for the data; default is 1.4826
    :type data: array-like
    :type nfact: float
    :returns nmad: (normalized) median absolute deviation of data.
    """
    m = np.nanmedian(data)
    return nfact * np.nanmedian(np.abs(data - m))

def nmad_outlier_removal(bins, DEM, dDEM, nfact=3):
    new_ddem = np.zeros(dDEM.shape)
    digitized = np.digitize(DEM, bins)
    for i, _ in enumerate(bins):
        if i == 0:
            print('Excluding outlier filtering on elevations lower than ' + str(bins[0]))
            this_bindata = dDEM[digitized == i]
            new_ddem[digitized == i] = this_bindata
            continue
        this_bindata = dDEM[digitized == i]
        this_nmad = nmad(this_bindata)
        this_bindata[np.abs(this_bindata) > nfact * this_nmad] = np.nan
        new_ddem[digitized == i] = this_bindata
    return new_ddem

# remove outliers 3-sigma filter
n_ddem = outlier_removal(bins, dem, ddem, nsig=3)
n_ddem_outside_glaciers = outlier_removal(bins, dem, ddem_outside_glaciers)

sig_ddem_out = np.where(np.isnan(n_ddem), n_ddem_outside_glaciers, n_ddem)

# out dir and filename
bname = os.path.basename(fp_ddem)
outname = bname[:-4] + '_edit_3sig_filtered.tif'
fp_out = out_dir + outname

# update profile and save
with rio.Env():
    
    # Update source profile: band count 1, set the
    # dtype to float32, specify LZW compression
    p_ddem.update(
        dtype=rio.float32,
        count=1,
        compress='lzw')
    
    # write new ddem to file
    with rio.open(fp_out, 'w', **p_ddem) as dst:
        dst.write(sig_ddem_out.astype(rio.float32), 1)



# normalized median abdolute difference filter
nmad_ddem = nmad_outlier_removal(bins, dem, ddem)
nmad_ddem_outside_glaciers = nmad_outlier_removal(bins, dem, ddem_outside_glaciers)
# fill areas outside glaciers with the ddem values
ddem_out = np.where(np.isnan(nmad_ddem), nmad_ddem_outside_glaciers, nmad_ddem)

#plt.imshow(nmad_ddem, cmap='BuPu')
#plt.show()

# 3 sigma thresholds
#thr_up = mean + stdev * 3
#thr_dn = mean - stdev * 3

# boolean mask 
#ddem_bool = np.logical_and(ddem < thr_up, ddem > thr_dn)
# select based on bool mask
#new_ddem = np.where(ddem_bool == True, ddem, np.nan)

# out dir and filename
bname = os.path.basename(fp_ddem)
outname = bname[:-4] + '_nmad_filtered.tif'
fp_out = out_dir + outname

# update profile and save
with rio.Env():
    
    # Update source profile: band count 1, set the
    # dtype to float32, specify LZW compression
    p_ddem.update(
        dtype=rio.float32,
        count=1,
        compress='lzw')
    
    # write new ddem to file
    with rio.open(fp_out, 'w', **p_ddem) as dst:
        dst.write(ddem_out.astype(rio.float32), 1)










