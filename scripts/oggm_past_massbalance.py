#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 07:43:00 2021

Use OGGM to get past mass balance. 

Input: edited glacier outlines and DEM for the starting year

Based on OGGM tutorials 'Getting started' and 'Run OGGM with GCM data'

Notes:
    Note that OGGM uses hydrological year 
    
    13.2.2021
    OGGM v1.4 latest build works on python 3.7 (MacOs)
    
    13.2.2021
    workflow.download_ref_tstars(base_url=params_url) <- workflow.py required editing,
    Added parameter follow_symlinks=True 
    
    13.2.2021
    Datasets in WGS84 (EPSG:4326)
    
To Do:
    
    
@author: apj
"""

import os
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import oggm
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import gis
from oggm.core import gcm_climate
from pyproj import CRS
from functools import partial, wraps
from shapely.ops import transform as shp_trafo
from salem import wgs84
from rgitools.funcs import compute_intersects
import shutil

# fig fp
figfp = r'/Users/apj/Documents/_HY/Greenland/OGGM'

# output path for fitted MB dataframe
mb_out = os.path.join(figfp, 'mb_result_fitted.csv')

cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border'] = 10

# read dem
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/study_area_dems_coregistered/80s_dem_studyarea_wgs84.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1953_masked_2212_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+1.82_y-13.07_z-2.21_align.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/1953_dem_wgs84.tif'
custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/elev_rasters/2016_dem_wgs84_studyarea.tif'

# geodetic MB
fp_geod = r'/Users/apj/Documents/_HY/Greenland/dVol/Outlines_dVolGtmwe_global.shp'
geodmb = gpd.read_file(fp_geod)

custom_dem_path

# set custom dem path to PATHS, so it can be used in tasks. tasks.define_glacier_region parameter source='user' will read dem from cfg.PATHS['dem_file']
cfg.PATHS['dem_file'] = custom_dem_path

# some parameters
cfg.PARAMS['use_intersects'] = True
cfg.PARAMS['use_rgi_area'] = False
cfg.PARAMS['continue_on_error'] = True

# read shapefile
fp = r'/Users/apj/Documents/_HY/Greenland/OGGM/outlines/2016_outlines_wgs84_singlepart.shp'
#fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/1953_extent_utm_edited_final_edit2_wgs84.shp'

# read outlines 
entity = gpd.read_file(fp)

# exclude glaciers smaller than 0.1 km2
entity = entity[entity['Area'] >= 1]

# geodetic MB list
geod_list = list(geodmb.RGIId)
# glaciers that have raised errors
rm_list = ['RGI60-05.02066', 'RGI60-05.02122', 'RGI60-05.02215', 'RGI60-05.02303', 'RGI60-05.01974', 'RGI60-05.02006', 'RGI60-05.02018', 'RGI60-05.02019', 'RGI60-05.02034', 'RGI60-05.02048', 'RGI60-05.02113' 'RGI60-05.02125', 'RGI60-05.02205', 'RGI60-05.02228', 'RGI60-05.02229', 'RGI60-05.02240', 'RGI60-05.02241', 'RGI60-05.02262', 'RGI60-05.02274', 'RGI60-05.02304', 'RGI60-05.02320', 'RGI60-05.02273_2', 'RGI60-05.02328_2', 'RGI60-05.01987_2']
rm_list2 = ['RGI60-05.01979','RGI60-05.02126_2', 'RGI60-05.01996', 'RGI60-05.02312_2', 'RGI60-05.02299','RGI60-05.01974', 'RGI60-05.02006','RGI60-05.02018', 'RGI60-05.02019', 'RGI60-05.02034', 'RGI60-05.02048', 'RGI60-05.02103', 'RGI60-05.02113', 'RGI60-05.02122', 'RGI60-05.02313', 'RGI60-05.02320','RGI60-05.02310_2', 'RGI60-05.02328_2', 'RGI60-05.01987_2', 'RGI60-05.02135']

# select only glaciers with geodetic MB
entity = entity[entity['RGIId'].isin(geod_list)]

# exclude
entity = entity[~entity.RGIId.isin(rm_list)]
entity = entity[~entity.RGIId.isin(rm_list2)]

# test with one glacier
#entity = entity[entity.RGIId == 'RGI60-05.01920']


# compute intersects
new_intersects = compute_intersects(entity)
# store intersects to working dir
entity_intersects_path = os.path.join(cfg.PATHS['working_dir'], 'entity_intersects.shp')
new_intersects.to_file(entity_intersects_path)

# set intersect file to use
cfg.set_intersects_db(new_intersects)

# year where outlines are valid
outline_year = os.path.basename(fp)[:4]

# check geometries
#for idx, row in entity.iterrows():
#    if entity.geometry.iloc[idx].type != 'Polygon':
#        print(row['RGIId'] + row.geometry.type)

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('User_past_mb', reset=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

#cfg.PATHS['working_dir'] = utils.gettempdir('user', reset=True)
#cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-GettingStarted', reset=True)

# glacier directory
gdirs = workflow.init_glacier_directories(entity)#, prepro_base_url=base_url) # <- from prepro_level uses preprocessed files
#rgi_id = gdir.rgi_id
#gdir.rgi_date = outline_year

# update gdir outline year
for gdir in gdirs:
    gdir.rgi_date = outline_year

# glacier region
#tasks.define_glacier_region(gdir, source='USER') # set user dem, this works for single glacier
workflow.execute_entity_task(tasks.define_glacier_region, gdirs, source='USER') # when using multiple glaciers

# plot glacier
#graphics.plot_domain(gdirs[0]) 

# print gdir files
#print(os.listdir(gdirs[0].dir))

# create glacier mask
#workflow.execute_entity_task(tasks.glacier_masks, gdirs)

# tasks to be executed
list_tasks = [
         tasks.glacier_masks, 
         tasks.compute_centerlines,
         tasks.initialize_flowlines,
         tasks.catchment_area,
         tasks.catchment_width_geom,
         tasks.catchment_width_correction,
         tasks.compute_downstream_line,
         tasks.compute_downstream_bedshape
         ]
for task in list_tasks:
    # some glaciers might result in error in some of the tasks, so use try to test if an exception is raised
    
    # The order matters!
    workflow.execute_entity_task(task, gdirs)
    # if the gdir raises exception, continue. The gdir raising the error should be removed

# list gdirs that raised error (manual task)
#rm_list = ['RGI60-05.02066', 'RGI60-05.02122', 'RGI60-05.02215', 'RGI60-05.02303']
# remove gdirs that raised error
#for gdir in gdirs:
#    if gdir.rgi_id in rm_list:
#        gdirs.remove(gdir)

# plot catchment areas
#graphics.plot_catchment_areas(gdirs, figsize=(8, 7))

### Climate tasks ###
# get tstars t* list and associated model parameters
params_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'

workflow.download_ref_tstars(base_url=params_url)

# run climate related entity tasks
workflow.climate_tasks(gdirs) # Downloads some files on the first time!

# list gdirs that raised error (manual task)
#rm_list = ['RGI60-05.01974', 'RGI60-05.02006', 'RGI60-05.02018', 'RGI60-05.02019', 'RGI60-05.02034', 'RGI60-05.02048', 'RGI60-05.02113' 'RGI60-05.02125', 'RGI60-05.02205', 'RGI60-05.02228', 'RGI60-05.02229', 'RGI60-05.02240', 'RGI60-05.02241', 'RGI60-05.02262', 'RGI60-05.02274', 'RGI60-05.02304', 'RGI60-05.02320', 'RGI60-05.02273_2', 'RGI60-05.02328_2', 'RGI60-05.01987_2']
# remove gdirs that raised error
#for gdir in gdirs:
#    if gdir.rgi_id in rm_list:
#        gdirs.remove(gdir)

# define year range
years = np.arange(1903, 2020)

# create dataframe to store results
mb_result = pd.DataFrame()

# Flowline Mass Balance
from oggm.core.massbalance import MultipleFlowlineMassBalance, PastMassBalance
for gdir in gdirs:
    rgi_id = gdir.rgi_id
    mbmod = MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True, mb_model_class=PastMassBalance)
    mb_ts = mbmod.get_specific_mb(year=years) # get mass balance
    # create dataframe of mb per year
    temp_df = pd.DataFrame({'mb':mb_ts}, index=years)
   
    # read geodetic mb for the glacier    
    test = geodmb[geodmb['RGIId'] == rgi_id]

    # get mm w.e. per year
    mmwe5385 = test['dmwe_53_85'].loc[test.index[0]] * 1000 / 32
    mmwe8516 = test['dmwe_85_16'].loc[test.index[0]] * 1000 / 31

    # avg difference to oggm mb per period. # Note hydrological year in OGGM. 
    davg5385 = np.average(temp_df['mb'].loc[1954:1986]) - mmwe5385
    davg8516 = np.average(temp_df['mb'].loc[1986:2017]) - mmwe8516
    
    # add difference and compute corrected mb to columns
#    temp_df['d5385'] = davg5385
#    temp_df['d8516'] = davg8516
    temp_df['corr_avg5385'] = temp_df['mb'].loc[1954:1985] - davg5385
    temp_df['corr_avg8516'] = temp_df['mb'].loc[1986:2017] - davg8516
    
    # join columns and drop the ones not needed
    temp_df['corr_mb'] = temp_df['corr_avg5385'].fillna(temp_df['corr_avg8516'])
    temp_df = temp_df.drop(['corr_avg5385', 'corr_avg8516'], axis=1)
    
    # compute bias (i.e. difference to geodetic mean) to column
    temp_df['bias'] = temp_df.mb - temp_df.corr_mb
    # add id
    temp_df['RGIId'] = gdir.rgi_id
    
    # append temporary df 
    mb_result = mb_result.append(temp_df)    

# save mb_results before grouping
mb_result.to_csv(mb_out, sep=';')

# group by year and average
mb_grouped = mb_result.groupby(mb_result.index).mean()

# calculate rolling average to new column
mb_grouped['avg_10year'] = mb_grouped['mb'].rolling(10,min_periods=1).mean()
mb_grouped['corr_avg_10year'] = mb_grouped['corr_mb'].rolling(10,min_periods=1).mean()
mb_grouped['avg_30year'] = mb_grouped['mb'].rolling(30,min_periods=1).mean()

# save dataframe
outname = os.path.join(figfp, 'OGGM_2016_calculated_MB_1902_2019_over1km2_GeodMB_fit_11052021.csv')
mb_grouped.to_csv(outname, sep=';')

# read saved dataframe csv
#mb_grouped = pd.read_csv(outname, sep=';')
# set index
#mb_grouped = mb_grouped.set_index('Unnamed: 0')

# compute and print averages
print('SMB in 1953-1985 ' + str(np.nanmean(mb_grouped['mb'].loc[1954:1986])))
print('SMB in 1985-2016 ' + str(np.nanmean(mb_grouped['mb'].loc[1986:2017])))
print('SMB in 1953-2016 ' + str(np.nanmean(mb_grouped['mb'].loc[1954:2017])))

print('SMB fit to geod in 1953-1985 ' + str(np.nanmean(mb_grouped['corr_mb'].loc[1954:1986])))
print('SMB fit to geod in 1985-2016 ' + str(np.nanmean(mb_grouped['corr_mb'].loc[1986:2017])))
print('SMB fit to geod in 1953-2016 ' + str(np.nanmean(mb_grouped['corr_mb'].loc[1954:2017])))

print('SMB in 1958-1996 ' + str(np.nanmean(mb_grouped['mb'].loc[1959:1997])))
print('SMB in 1997-2015 ' + str(np.nanmean(mb_grouped['mb'].loc[1998:2016])))


# plot
fig, ax = plt.subplots()
ax.plot(mb_grouped.index, mb_grouped.mb, color='#CA6F1E', alpha=0.7, linewidth=1, label='OGGM MB')
ax.plot(mb_grouped.index, mb_grouped.corr_mb, color='#249CF4', linewidth=1, label='Fitted to geodetic MB')
ax.plot(mb_grouped.index, mb_grouped.avg_10year, color='k', ls='-.', linewidth=1.5, alpha=0.5, label='10-year average')
ax.plot(mb_grouped.index, mb_grouped.corr_avg_10year, color='k', ls='-.', linewidth=1.5, label='Fitted 10-year average')

#ax.plot(mb_grouped.index, mb_grouped.avg_30year, color='g', label='30-year average')
ax.axvline(1954, ls='--', color='k', alpha=0.5)
ax.axvline(1986, ls='--', color='k', alpha=0.5)
ax.axvline(2017, ls='--', color='k', alpha=0.5)
ax.set_ylabel('SMB (mm $yr^-$$^1$)')
ax.legend()
fig.suptitle('OGGM mass balance 1902-2019')
fig.tight_layout()

plt.savefig(os.path.join(figfp, 'OGGM_calculated_MB_1902-2019_geodMBFit.png'), dpi=150)











