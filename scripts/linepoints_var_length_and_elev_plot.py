#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 14:26:34 2020

Plotting glacier centerline elevations from multiple line point features with elevation attribute

Notes:
    21.10.2020
    Multiple inputs where to plot centerlines
    
    
    Tributary glaciers should have unique number to prevent overwriting plots

@author: apj
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import os
from shapely.geometry import Point, LineString
from scipy.spatial import distance
import numpy as np
import ogr
import glob


# set filepath for folder with shapefiles
fp = r'/Users/apj/Documents/_HY/Greenland/centerlines/edited_glacier_divides/points_w_elev/*.shp'


# function to get coordinate list
def getCoords(dframe):
    # empty list for coordinates
    coord_list = []
    for i in dframe['geometry']:
        p = str(i)
        #print(p)
        point = ogr.CreateGeometryFromWkt(p)
        xy = point.GetX(), point.GetY()
        coord_list.append(xy)
    return coord_list

# function to calculate line length
def calculateLength(coordinate_list):
    # calculate distances between points to a list
    d = 0
    lengths = []
    while d < len(coordinate_list):
        if d == 0:
            dist = 0
        elif d >= 1:
            dist = distance.euclidean(coordinate_list[d-1], coordinate_list[d])
        elif d == len(coordinate_list):
            break
        #print(dist)
        lengths.append(dist)        
        d += 1
    # return list of euclidean distances between points
    return lengths



# empty list for files
flist = []

# read files to list
for fpath in glob.glob(fp):
    f = gpd.read_file(fpath)
    # check column 'MAIN' type and change to  float if not already
    f = f.astype({'MAIN': float})
    flist.append(f)    


# get unique glacier centerline IDs to list. Unique IDs should be same between input files.
# Could add a check that unique IDs match between input files
unique_ids = list(flist[0]['RGIID2'].unique())


"""
# combine unique RGI ids from dataframes
for df in flist:
    # new column with RGIID and ORIG_FID
    df['RGIID_Tr'] = df['RGIID'] + '_' + df['MAIN'].astype(str)

    # get unique field id's to list
    unique_orig_fids = list(df['RGIID_Tr'].unique())
        
    # add unique ids to list
    [all_unique_orig_fids.append(x) for x in unique_orig_fids if x not in all_unique_orig_fids] 
"""

# function to calculate line cumulative length to new column
def calcLength(selected_dataframe):
    # get coordinates
    coords = getCoords(selected_dataframe)
    
    # get length between points
    length_between_points = calculateLength(coords)
    
    # assign lengths to new column in dataframe subselection
    selected_dataframe = selected_dataframe.assign(lengths=length_between_points)
    
    # compute cumulative length
    selected_dataframe['cum_length'] = np.cumsum(length_between_points)
    # length in km
    selected_dataframe['cum_length_km'] = selected_dataframe['cum_length'] / 1000
    
    return selected_dataframe


# counter
counter = 1

# loop through unique IDs and create plots. NOTE!! check df order in list and elevation column
for uniq_id in unique_ids:
    
    # selection from the dataframes by Unique id
    sel1 = flist[0][flist[0]['RGIID2'] == uniq_id]
    sel2 = flist[1][flist[1]['RGIID2'] == uniq_id]
    sel3 = flist[2][flist[2]['RGIID2'] == uniq_id]
    
    df1 = calcLength(sel1)
    df2 = calcLength(sel2)
    df3 = calcLength(sel3)
       
    # plot elevations and glacier length with matplotlib
    plt.plot(df1['cum_length_km'], df1['2010s_dem'], 'k-')
    plt.plot(df2['cum_length_km'], df2['50s_dem'], 'b-')
    plt.plot(df3['cum_length_km'], df3['80s_dem'], 'r-')
    #plt.plot(line_subset['cum_length'], line_subset['ArcticDEM_'], 'g-')
    if df1['MAIN'].all() != 1:
        plt.title('Glacier centerline elevations for ' + str(df1['RGIID'].iloc[0]) + ' tributary glacier')
    else:
        plt.title('Glacier centerline elevations for glacier: ' + str(df1['RGIID'].iloc[0]) + ' main trunk')
    plt.xlabel('distance (km)')
    plt.ylabel('elevation (m)')
    plt.legend(['2014 DEM', '50s DEM', 'AERODEM'], loc='upper right')
    plt.show
        
    """
    Saving figure
    """
    # define filename and filepath for the figure and save figure
    figname = str(uniq_id) + '.png'
    out_path = r'/Users/apj/Documents/_HY/Greenland/centerlines/figures/edited_divides_w_tributaries/'
    figpath = os.path.join(out_path, figname)
    plt.savefig(figpath, format = 'png')
    
    # pause before closing the figure
    plt.pause(0.5)
    plt.close()
    print(str(counter) + ' / ' + str(len(unique_ids)))
    counter += 1















