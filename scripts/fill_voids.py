#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:36:47 2020

Fill voids using elevaiton bins in differenced DEM. 

Inputs:
    glacier outline shapefile
    DEM where to get elevation bins
    dDEM - differenced DEM

Elevation bins, mean elevation difference per elevation bin per glacier (local hypsometric)
Test conditions:
   - more than 40% of bin area has valid values - if not: reject bin (fill later with polynomial)
   - entire glacier has more than 2/3 of elevation range covered - if not: exclude glacier
   - less than 1/3 of glacier area has valid values: exclude glacier
        
Glacier ID of the glaciers to exclude are saved to a list. 

To Do:
    clean the script

@author: apj. 28.12.2020
"""

import fiona
import numpy as np
from numpy.polynomial.polynomial import polyval, polyfit
import rasterio as rio
import rasterio.mask
import pandas as pd
import geopandas as gpd
import glob
import matplotlib.pyplot as plt
import os

# filepath
fp_shape = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_final_edit.shp'
fp_tif = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filtered/*nmad_filtered.tif'
#fp_dem = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_22122020/1953_masked_2212_2016_diss_ext_masked_nuth_x-2.25_y-18.73_z-29.21_align.tif'
fp_dem2 = r'/Users/apj/Documents/_HY/Greenland/kfs_dem/2016dem_wgs84_ellipsoid_utm_25m_studyarea_edit.tif'
fp_gl_mask = r'/Users/apj/Documents/_HY/Greenland/masks/glacier_mask2.tif'

# select void fill method, local hypsometric or global
method = 'global'


# read file
gdf = gpd.read_file(fp_shape)
# list tiff files
tiffiles = []
for tiffp in glob.glob(fp_tif):
    tiffiles.append(tiffp)

# read dem and create elevation bins
#with rio.open(fp_dem) as src:
#    dem = src.read(1)
#    dem_nodata = src.nodata
    
# change nodata to np.nan
#dem[(dem == dem_nodata)] = np.nan

# read shapefile geometries
with fiona.open(fp_shape, "r") as shapefile:
    shapes = [feature["geometry"] for feature in shapefile]

# mask dem to find out glacier elevation range
#with rio.open(fp_dem) as src:
#    out_img, out_transform = rasterio.mask.mask(src, shapes, crop=False)
#    out_nodata = src.nodata

#elev_range = out_img[0]
#elev_range[(elev_range == out_nodata)] = np.nan

with rio.open(fp_dem2) as src:
    dem = src.read(1)
    dem_nodata = src.nodata

dem[dem == dem_nodata] = np.nan


# get elevation bins
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


# function to create glacier mask from non-voided raster and shapefile
def createGlacierMask(fp_dem2, shapefile):    
    # create glacier mask for ddem
    with rasterio.open(fp_dem2) as src:
        glac_mask, glac_out_transform = rasterio.mask.mask(src, shapes, crop=False)
        glac_out_meta = src.meta
        glac_nodata = src.nodata
        glac_p = src.profile
    
    # select array
    glac = glac_mask[0]
    glac[(glac == glac_nodata)] = np.nan
    glacier_mask = np.isfinite(glac)
    
    gl_mask_bin = np.where(glacier_mask == False, 0, 1)
    
    # write to file
    # update profile and save
    with rio.Env():
        
        # Update source profile: band count 1, set the
        # dtype to float32, specify LZW compression
        glac_p.update(
            dtype=rio.int8,
            nodata = 0,
            count=1,
            compress='lzw')
        
        # write new ddem to file
        with rio.open(fp_gl_mask, 'w', **glac_p) as dst:
            dst.write(gl_mask_bin.astype(rio.int8), 1)



# functions from pybob
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

if method == 'local':
    
    # empty list for glaciers to be excluded
    exclude_glaciers = []
    
    # ddem filepath
    fp_ddem = tiffiles[1]
    
    # new ddem to store results from single glaciers
    new_ddem = np.zeros(dem.shape)
    
    # empty array to store all results
    out_ddem = np.empty(dem.shape)
    out_ddem[:] = np.nan
    
    counter = 1
        
    # loop through glaciers, mask ddem and fill voids
    for glac_i in gdf['RGIId']:
        # select glacier
        sel = gdf[gdf['RGIId'] == glac_i] # test glacier 'RGI60-05.01916'
        gl_id = sel['RGIId'][sel['RGIId'].index[0]]
        #break
        
        # get ddem data within glacier
        with rasterio.open(fp_ddem) as src:
            #ddem = src.read(1)
            out_image, out_transform = rasterio.mask.mask(src, sel.geometry, crop=False) # mask with polygon
            ddem_nodata = src.nodata
            p_ddem = src.profile
            out_meta = src.meta
    
        # ddem of glacier area and set outlying area to nodata
        ddem = out_image[0]
        ddem[(ddem == ddem_nodata)] = np.nan
    
        # glacier mask
        with rasterio.open(fp_gl_mask) as src:
            #ddem = src.read(1)
            mask_image, mask_transform = rasterio.mask.mask(src, sel.geometry, crop=False) # mask with polygon
            mask_nodata = src.nodata
            p_mask = src.profile
            
        gl_mask = mask_image[0]
        gl_mask = np.where(gl_mask == 1, True, False)    
        
        # valid area mask
        valid = np.logical_and(gl_mask, np.logical_and(np.isfinite(dem), np.isfinite(ddem)))
        whole = gl_mask[gl_mask == True] # where glacier mask is True
            
        # extract valid area
        dem_data = dem[valid]
        ddem_data = ddem[valid]
        
        # check that array is not empty
        if ddem_data.size == 0:
            print('Zero sized array on glacier: ' + glac_i)
            counter += 1
            continue
        
        # get elevation bins for the selected glacier 
        gl_bins = getBins(dem_data, 100)
        
        # classify dem to bins
        digitized = np.digitize(dem, gl_bins)
        
        # Step 1.1
        # loop through bins, select ddem data in the bin and check data coverage in the bin
        for i, _ in enumerate(gl_bins):
               
            # get data by bin
            this_bindata = ddem[digitized == i]
            this_binmask = gl_mask[digitized == i]
            
            #print(gl_bins[i])
            #print(np.count_nonzero(~np.isnan(this_bindata)))
            #print(np.count_nonzero(this_binmask))
    
            # 1st threshold
            # check percentage of valid values in the bin inside the glacier outline
            n_values_this_bindata = np.count_nonzero(~np.isnan(this_bindata))
            n_values_this_binmask = np.count_nonzero(this_binmask)
            
            # if valid data in the bin less than  %, exclude bin
            if n_values_this_bindata < ( 0.4 * n_values_this_binmask):
                print('Rejecting bin: ' + str(gl_bins[i]))
                this_bindata = np.full_like(this_bindata, np.nan)
        
        # place data to bins and calculate statistics
        bin_mean = bin_data(gl_bins, ddem_data, dem_data, mode='mean', nbinned=False)
       
        # 2nd threshold
        # check that the data in elevation bins covers more than 2/3 of the glacier elevation range
        _bins = gl_bins[np.isfinite(bin_mean)] # get valid bins
        if len(_bins) < (2/3 * len(gl_bins)):
            # add to excluded list
            print('Less than 2/3 coverage in binned data on glacier: ' + gl_id)
            exclude_glaciers.append(gl_id)
        
        # polyfit
        def fitPoly(bin_means):
            _bins = gl_bins[np.isfinite(bin_means)] # get valid bins
            _curve = bin_means[np.isfinite(bin_means)] # get mean values from valid bins
            p = polyfit(_bins, _curve, 3) # fit polynomial
            fill_ = polyval(gl_bins, p) # polynomial value
            bin_means[np.isnan(bin_means)] = fill_[np.isnan(bin_means)] # fill nan bin
            return bin_means 
           
        #### Step 2 - check valid area covered by ddem_data
        
        # compare valid data and whole data
        #prop_covered = ddem_data.size / whole.size
        
        #if prop_covered < (1/3 * whole.size):
        if np.count_nonzero(ddem_data) < (1/3 * np.count_nonzero(gl_mask)):
            print('Less than one third coverage on glacier: ' + gl_id)
            exclude_glaciers.append(gl_id)
        
        # Step 3
    #    temp_bin_mean = np.full_like(this_bindata, bin_means[i])
        # fill void cells in the bins
    #    filled_ddem = np.where(np.logical_and(np.isnan(this_bindata), this_binmask==True), temp_bin_mean, this_bindata)
    
        filled_ddem = np.zeros(ddem.shape)
        # fill values in the bins
        for i, _ in enumerate(gl_bins):
            # if bin is empty, fit polynomial and get value
            if np.isnan(bin_mean[i]):
                bin_poly = fitPoly(bin_mean)
                filled_ddem[digitized == i] = bin_poly[i]
                print('Polyfill bin ' + str(_))
            else: # assign bin mean to voids   
                # get data by bin
                test_bindata = filled_ddem[digitized == i]
                this_bindata_ddem = ddem[digitized == i] 
                
                # filled bindata
                filled_bindata = np.where(np.isnan(this_bindata_ddem), bin_mean[i], test_bindata)
        
                filled_ddem[digitized == i] = filled_bindata        
                
        # update ddem
        fill_area = np.logical_and(gl_mask, np.isnan(ddem))
        
        # filled glacier
        new_ddem = np.where(fill_area == True, filled_ddem, ddem)
        
        # add single glacier result in new_ddem to out_ddem
        out_ddem = np.where(np.isfinite(out_ddem), out_ddem, new_ddem)
        
        #plt.imshow(out_ddem)
        #plt.show()
        
        # print progress
        print('Done with ' +  str(counter) + ' / ' + str(len(gdf)))
        
        #if counter == 10:
        #    break
        counter += 1
            
    # fill areas outside glacier with ddem values (optional)
    #with rio.open(fp_ddem) as src:
    #    ddem_fill = src.read(1)
    
    #out_ddem = np.where(np.isfinite(out_ddem), out_ddem, ddem_fill)

    # output filename and write filled ddem to file
    fill_dem_dir = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filled_ddem/'
    bname = os.path.basename(fp_ddem)
    fname = bname[0:12] + '_local_hyps_filled_ddem.tif'
    filled_ddem_out = fill_dem_dir + fname

if method == 'global':
    
    for tif in tiffiles:
        
        fp_ddem = tif
         # new ddem to store results 
        new_ddem = np.zeros(dem.shape)
        # empty array for bin means
        filled_ddem = np.zeros(dem.shape)
        
        # empty array to store all results
        out_ddem = np.empty(dem.shape)
        out_ddem[:] = np.nan
        
        counter = 1
        
         # get ddem data within shapes
        with rasterio.open(fp_ddem) as src:
            #ddem = src.read(1)
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=False) # mask with polygon
            ddem_nodata = src.nodata
            p_ddem = src.profile
            out_meta = src.meta
    
        # ddem of glacier area and set outlying area to nodata
        ddem = out_image[0]
        ddem[(ddem == ddem_nodata)] = np.nan
        
        # glacier mask
        with rasterio.open(fp_gl_mask) as src:
            #ddem = src.read(1)
            mask_image, mask_transform = rasterio.mask.mask(src, shapes, crop=False) # mask with polygon
            mask_nodata = src.nodata
            p_mask = src.profile
            
        gl_mask = mask_image[0]
        gl_mask = np.where(gl_mask == 1, True, False)    
        
        # valid area mask
        valid = np.logical_and(gl_mask, np.logical_and(np.isfinite(dem), np.isfinite(ddem)))
        whole = gl_mask[gl_mask == True] # where glacier mask is True
            
        # extract valid area
        dem_data = dem[valid]
        ddem_data = ddem[valid]
        
        # check that array is not empty
        if ddem_data.size == 0:
            print('Zero sized array on glacier: ' + glac_i)
            counter += 1
            continue
        
        # get elevation bins for the selected glacier 
        gl_bins = getBins(dem_data, 100)
        
        # classify dem to bins
        digitized = np.digitize(dem, gl_bins)
        
        # place data to bins and calculate statistics
        bin_mean = bin_data(gl_bins, ddem_data, dem_data, mode='mean', nbinned=False)
       
        # fill values in the bins
        for i, _ in enumerate(gl_bins):
            # if bin is empty, fit polynomial and get value
            if np.isnan(bin_mean[i]):
                print('Empty bin')
                #bin_poly = fitPoly(bin_mean)
                #filled_ddem[digitized == i] = bin_poly[i]
                #print('Polyfill bin ' + str(_))
            else: # assign bin mean to voids   
                # get data by bin
                test_bindata = filled_ddem[digitized == i]
                this_bindata_ddem = ddem[digitized == i] 
                
                # filled bindata
                filled_bindata = np.where(np.isnan(this_bindata_ddem), bin_mean[i], test_bindata)
        
                filled_ddem[digitized == i] = filled_bindata        
                
        # update ddem
        fill_area = np.logical_and(gl_mask, np.isnan(ddem))
        
        # filled glacier
        new_ddem = np.where(fill_area == True, filled_ddem, ddem)
        
        # output filename and write filled ddem to file
        fill_dem_dir = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filled_ddem/'
        bname = os.path.basename(fp_ddem)
        fname = bname[0:12] + '_global_hyps_filled_ddem.tif'
        filled_ddem_out = fill_dem_dir + fname
        
        
        # update profile and save
        with rio.Env():
            
            # Update source profile: band count 1, set the
            # dtype to float32, specify LZW compression
            p_ddem.update(
                dtype=rio.float64,
                count=1,
                compress='lzw')
            
            # write new ddem to file
            with rio.open(filled_ddem_out, 'w', **p_ddem) as dst:
                dst.write(new_ddem.astype(rio.float64), 1)


        
# save exclude glacier list to file
#out_lis = r'/Users/apj/Documents/_HY/Greenland/dem_diff/filled_ddem_bin_thresh05/glaciers_to_exclude.csv'
gdf_out = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_exclude_bin_thresh.shp'
#excluded_series = pd.Series(exclude_glaciers)
#excluded_series.to_csv(out_list, name='IDs')
gdf_update = gdf[~gdf['RGIId'].isin(exclude_glaciers)]
gdf_update.to_file(gdf_out, driver='ESRI Shapefile')

# update profile and save
with rio.Env():
    
    # Update source profile: band count 1, set the
    # dtype to float32, specify LZW compression
    p_ddem.update(
        dtype=rio.float64,
        count=1,
        compress='lzw')
    
    # write new ddem to file
    with rio.open(filled_ddem_out, 'w', **p_ddem) as dst:
        dst.write(out_ddem.astype(rio.float64), 1)








