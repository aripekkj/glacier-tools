#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 18:35:22 2020

Calculate volume differences and uncertainties


@author: apj
"""

import geopandas as gpd
import pandas as pd
import math
from rasterstats import zonal_stats
import fiona
import glob
import os
import numpy as np
import rasterio as rio
import rasterio.mask
import matplotlib.pyplot as plt

fp_glaciers = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_final_edit.shp'
#fp_glaciers_out = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/Dissolved_outline_50s80s2010s_utm_final_edit_excluded.shp'
fp_tiff = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filled_ddem/*linear*.tif'
fp_filtered = r'/Users/apj/Documents/_HY/Greenland/dem_diff/vgridshift/filtered/*nmad*.tif'
fp_out = r'/Users/apj/Documents/_HY/Greenland/dVol/Outlines_dVol_bin_thresh.shp'
fp_exclude = r'/Users/apj/Documents/_HY/Greenland/dem_diff/filled_ddem/glaciers_to_exclude_edit.csv'
# co-register surfaces
fp_53surface = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1953_masked_2212_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+1.82_y-13.07_z-2.21_align_filt.tif'
fp_85surface = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1985_rm401_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+24.98_y-4.86_z-1.84_align_filt.tif'

def excludeByID(excluded_list, in_gdf, id_column):
    # select all but ones in the list
    selection = in_gdf[~in_gdf[id_column].isin(excluded_list)]
    return selection

# read shapefile geometries
def readGeom(fp_shapefile):    
    with fiona.open(fp_shapefile, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
        return shapes

# cell resolution
r = 25

# read file
glaciers = gpd.read_file(fp_glaciers)
#exclude = pd.read_csv(fp_exclude) # glaciers to exclude
# rename columns
#exclude = exclude.rename(columns={'Unnamed: 0': 'index', '0': 'ID'})
#remoce duplicates and drop extra columns
#exclude = exclude.drop_duplicates(subset='ID')
#exclude = exclude.drop(columns=['index'])
# to list
#excludelist = exclude['ID'].tolist()

# exclude and select certain glaciers 
#glaciers = excludeByID(excludelist, glaciers, 'RGIId')

# list tiff files
tiffiles = []
for file in glob.glob(fp_tiff):
    tiffiles.append(file)

filtered = []
for filtfile in glob.glob(fp_filtered):
    filtered.append(filtfile)

# root mean square error for a column of values
def getRMS(df, column_name):
    """
    
    Parameters
    ----------
    df : DataFrame
    column_name : Column name, String

    Returns
    -------
    rms : float
        Root mean square

    """
    avg = df[column_name].mean()
    rms = math.sqrt(avg**2)
    return rms

def ddemRMS(ddem, coreg_surface):
    """
    
    Parameters
    ----------
    fp_ddem : TYPE String
        DESCRIPTION. Raster filepath

    Returns
    -------
    RMS

    """        
    # read files
    with rio.open(ddem) as src:
        ddem = src.read(1)
        ddem_nodata = src.nodata
        
    with rio.open(coreg_surface) as src:
        surf = src.read(1)
        surf_nodata = src.nodata
        
    # make sure nodata is nan
    ddem[(ddem == ddem_nodata)] = np.nan
    surf[(surf == surf_nodata)] = np.nan
    
    # select cells from ddem
    ddem_sel = np.where(np.isfinite(surf), ddem, np.nan)
    
    # compute RMS for co-registered cells
    squared = ddem_sel**2
    mean_sq = np.nanmean(squared)
    rms_ddem = math.sqrt(mean_sq)
    return rms_ddem    

# rms
rms_53_2016 = ddemRMS(filtered[1], fp_53surface)
rms_85_2016 = ddemRMS(filtered[2], fp_85surface)
rms_53_85 = rms_53_2016 + rms_85_2016

# compute volume difference and uncertainties for whole area
glac = readGeom(fp_glaciers)
for tif in tiffiles:
    tif_bname = os.path.basename(tif)
    # read raster and mask with selected polygon
    with rio.open(tif) as src:
        out_image, out_transform = rasterio.mask.mask(src, glac, crop=True)
        img_nan = src.nodata

    ddem_data = out_image[0]
    ddem_data[(ddem_data == img_nan)] = np.nan
    
    meandh = np.nanmean(ddem_data)
    diffVol = np.nansum(ddem_data) * r**2
    diffVol_km3 = diffVol / 10**9

    #diffGt = diffVol_km3 * ice_density

    #print('Volume change in ' + tif_bname[0:12] + ': ' + str(round(diffGt, 3)) + ' Gt')
    print('Average elevation change in ' + tif_bname[0:12] + ': ' + str(round(meandh, 3)) + ' meters')
    
counter = 0
# compute error in volume per glacier
for i in glaciers['RGIId']:
    #print(i)
    
    # select one glacier
    sel = glaciers[glaciers['RGIId'] == i]
    # glacier area
    temp = sel.geometry
    area = temp.geometry.area[temp.index[0]]
    #area = row.geometry.area
    # compute mean elevation change 
    for tif in tiffiles:
        tif_bname = os.path.basename(tif)
        # read raster and mask with selected polygon
        with rio.open(tif) as src:
            out_image, out_transform = rasterio.mask.mask(src, temp.geometry, crop=True)
            img_nan = src.nodata
            p = src.profile
        
        # masked tif
        masked = out_image[0]
        
        # change no data to nan
        masked[(masked == img_nan)] = np.nan
        
        # count cells
        #cells = np.count_nonzero(masked) # probably wrong
        #cells = area / r / r # how many cells fits inside the area
        
        # plot
        #plt.imshow(masked, cmap='GnBu')
        #plt.show()
        
        # elevation change
        dHsum = np.nansum(masked)
        
        # calculate volume change (sum of elevation difference times xy resolution)
        dVol = np.nansum(masked) * r**2
        #dVolkm3 = dVol / 10**9 # convert to cubic kilometers
        
        # volume change in gigatonnes
        #dGt = dVol * ice_density
        
        # get mean elevation difference        
        mean_dh = np.nanmean(masked)
    
        # select which rms error to use
        if '1953' and '1985' in tif_bname:
            rand = rms_53_85
        elif '1953' and '2016' in tif_bname:
            rand = rms_53_2016
        else:
            rand = rms_85_2016
    
        
        # combine column name
        eVol_colname = 'eV_' + tif_bname[2:4] + '_' + tif_bname[10:12]
        dVol_colname = 'dV_' + tif_bname[2:4] + '_' + tif_bname[10:12]
        eGt_colname = 'eGt_' + tif_bname[2:4] + '_' + tif_bname[10:12]
        dGt_colname = 'dGt_' + tif_bname[2:4] + '_' + tif_bname[10:12]
        
        # store error to new column in shapefile
        #glaciers.loc[glaciers['RGIId'] == i, eVol_colname] = eVol
        glaciers.loc[glaciers['RGIId'] == i, dVol_colname] = dVol
        #glaciers.loc[glaciers['RGIId'] == i, eGt_colname] = eGt
        #glaciers.loc[glaciers['RGIId'] == i, dGt_colname] = dGt
        glaciers.loc[glaciers['RGIId'] == i, 'area_m2'] = area
    
    counter += 1
    #print(str(counter) + ' / ' + str(len(glaciers)))


# some constant values
L = 500 # autocorrelation distance
r = 25 # pixel size
eA = 0.1 # error in glacier area (percentage)

# ice density (Gt/km3)
ice_density = 0.900

# sum volume and area
Vsum = sum(glaciers['dV_53_16'])
Vsum = sum(glaciers['dV_53_85'])
Vsum = sum(glaciers['dV_85_16'])
Asum = sum(glaciers['area_m2'])

def getUncertainty(df, volume_column, area_column, rand_e):
    Volsum = sum(df[volume_column])
    Areasum = sum(df[area_column])
    
    cells = Asum / r / r # how many cells fits inside the area
    
    # elevation change
    dh = Volsum / Areasum
    
    # eqn. 2 from the paper - volume change uncertainty due to dh uncertainty
    e_h = np.sqrt(rand_e**2 / (np.sqrt(cells/(L/r)**2))) * Areasum
    # volume change uncertainty due to area uncertainty
    e_a = Volsum * eA
    
    # eqn. 3 from the paper - total uncertainty in volume change, per year
    total_uncert = np.sqrt(e_h**2 + e_a**2) / Areasum # total 
    print('Total volume change for ' + volume_column + ': {:.2f}±{:.2f} m'.format(dh, total_uncert))

    return dh, total_uncert

result_53_16 = getUncertainty(glaciers, 'dV_53_16', 'area_m2', rms_53_2016)
result_53_85 = getUncertainty(glaciers, 'dV_53_85', 'area_m2', rms_53_85)
result_85_16 = getUncertainty(glaciers, 'dV_85_16', 'area_m2', rms_85_2016)



# write result to file
#glaciers.to_file(fp_out, driver='ESRI Shapefile')









