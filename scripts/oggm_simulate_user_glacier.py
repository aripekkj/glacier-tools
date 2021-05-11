#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 07:43:00 2021

Run OGGM simulation. Based on OGGM tutorials 'Getting started' and 'Run OGGM with GCM data'

Notes:
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
figfp = r'/Users/apj/Documents/_HY/Greenland/OGGM/gcm_figs'

cfg.initialize(logging_level='WARNING')
cfg.PARAMS['border'] = 10

# read dem
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/study_area_dems_coregistered/80s_dem_studyarea_wgs84.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1953_masked_2212_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+1.82_y-13.07_z-2.21_align.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/1953_dem_wgs84.tif'
custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/2016_dem_wgs84_studyarea.tif'

custom_dem_path

# set custom dem path to PATHS, so it can be used in tasks. tasks.define_glacier_region parameter source='user' will read dem from cfg.PATHS['dem_file']
cfg.PATHS['dem_file'] = custom_dem_path

# some parameters
cfg.PARAMS['use_intersects'] = True
cfg.PARAMS['use_rgi_area'] = False
cfg.PARAMS['continue_on_error'] = True

# read shapefile
#fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/1953_extent_utm_edited_final_edit_wgs.shp'
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/2016_extent_updated_divides_utm_21102020_fixed_final_edit_wgs84.shp'

# list of glaciers to use/exclude
sel = ['RGI60-05.01979', 'RGI60-05.02299','RGI60-05.01974', 'RGI60-05.02006','RGI60-05.02018', 'RGI60-05.02019', 'RGI60-05.02034', 'RGI60-05.02048', 'RGI60-05.02103', 'RGI60-05.02113', 'RGI60-05.02122', 'RGI60-05.02313', 'RGI60-05.02320','RGI60-05.02310_2', 'RGI60-05.02328_2', 'RGI60-05.01987_2', 'RGI60-05.02135']

# read outlines 
entity = gpd.read_file(fp)

# exclude glaciers smaller than 0.1 km2
entity = entity[entity['Area'] >= 1]

# compute intersects
#new_intersects = compute_intersects(entity)
# store intersects to working dir
#entity_intersects_path = os.path.join(cfg.PATHS['working_dir'], 'entity_intersects.shp')
#new_intersects.to_file(entity_intersects_path)

# set intersect file to use
#cfg.set_intersects_db(new_intersects)

# select subset
entity = entity[~entity.RGIId.isin(sel)]

# year where outlines are valid
outline_year = os.path.basename(fp)[:4]

# check geometries
#for idx, row in entity.iterrows():
#    if entity.geometry.iloc[idx].type != 'Polygon':
#        print(row['RGIId'] + row.geometry.type)

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('User_gcm_data', reset=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

#cfg.PATHS['working_dir'] = utils.gettempdir('user', reset=True)
#cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-GettingStarted', reset=True)

# test with one glacier
#entity = entity[entity['RGIId'] == 'RGI60-05.02114']

# rgi_ids for default data
#rgi_ids = ['RGI60-05.01916']

# glacier directory
#base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/L3-L5_files/CRU/centerlines/qc3/pcp2.5/no_match/'
gdirs = workflow.init_glacier_directories(entity)#, prepro_base_url=base_url) # <- from prepro_level uses preprocessed files
#gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3, prepro_border=80) # this is loading default dataset
#rgi_id = gdir.rgi_id
#gdir.rgi_date = outline_year


#for gdir in gdirs:

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

        
# plot catchment areas
#graphics.plot_catchment_areas(gdirs, figsize=(8, 7))

### Climate tasks ###
# get tstars data to working dir
tstar_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'

workflow.download_ref_tstars(base_url=tstar_url)

# run climate related entity tasks
workflow.climate_tasks(gdirs) # Downloads some files on the first time!

# remove glacier that caused error, setting rgi IDs is manual
#for gdir in gdirs:
#    if gdir.rgi_id == 'RGI60-05.01510':
#        gdirs.remove(gdir)

# Flowline Mass Balance
from oggm.core.massbalance import MultipleFlowlineMassBalance
for gdir in gdirs:
    mbmod = MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)

# Ice thickness
list_talks = [
         tasks.prepare_for_inversion,  # This is a preprocessing task
         tasks.mass_conservation_inversion,  # This does the actual job
         tasks.filter_inversion_output  # This smoothes the thicknesses at the tongue a little
         ]
for task in list_talks:
    workflow.execute_entity_task(task, gdirs)

# Convert the flowlines to a "glacier" for the ice dynamics module
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs);

###  GCM data simulation ### 
# climate data filepaths
bp = 'https://cluster.klima.uni-bremen.de/~oggm/cmip5-ng/pr/pr_mon_CCSM4_{}_r1i1p1_g025.nc'
bt = 'https://cluster.klima.uni-bremen.de/~oggm/cmip5-ng/tas/tas_mon_CCSM4_{}_r1i1p1_g025.nc'
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
    # Download the files
    ft = utils.file_downloader(bt.format(rcp))
    fp = utils.file_downloader(bp.format(rcp))
    # bias correct them
    workflow.execute_entity_task(gcm_climate.process_cmip_data, gdirs, 
                                 filesuffix='_CCSM4_{}'.format(rcp),  # recognize the climate file for later
                                 fpath_temp=ft,  # temperature projections
                                 fpath_precip=fp,  # precip projections
                                 );

# run different scenarios
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
    rid = '_CCSM4_{}'.format(rcp)
    workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, ys=int(outline_year), 
                                 climate_filename='gcm_data',  # use gcm_data, not climate_historical
                                 climate_input_filesuffix=rid,  # use the chosen scenario
                                 output_filesuffix=rid,  # recognize the run for later
                                );

# DataFrames to store modelling results
rcp26_result = pd.DataFrame()
rcp45_result = pd.DataFrame()
rcp60_result = pd.DataFrame()
rcp85_result = pd.DataFrame()

# plot model results 
for gdir in gdirs:
    
    # plot modelling results
#    f, ax1 = plt.subplots(1, 1, figsize=(14, 4))
    for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
        rid = '_CCSM4_{}'.format(rcp)
        ds = utils.compile_run_output(gdir, input_filesuffix=rid)
   # store modelling results to dataframe
        temp_df = ds.to_dataframe()
        temp_df = temp_df[['calendar_year', 'volume', 'area', 'length']] # keep only relevant columns
        if rcp == 'rcp26':
            rcp26_result = rcp26_result.append(temp_df)
        elif rcp == 'rcp45':
            rcp45_result = rcp45_result.append(temp_df)
        elif rcp == 'rcp60':
            rcp60_result = rcp60_result.append(temp_df)
        elif rcp == 'rcp85':
            rcp85_result = rcp85_result.append(temp_df)

#        ds.isel(rgi_id=0).volume.plot(ax=ax1, label=rcp);
        #ds.isel(rgi_id=1).volume.plot(ax=ax2, label=rcp);
#    plt.legend();
#    plt.savefig(os.path.join(figfp, gdir.rgi_id + '_GCM_model_result_plot.png'), dpi=150)
    
    # Plot modelling output on map
#    f, (axs) = plt.subplots(4, 3, figsize=(14, 6))
    rn = 0
    for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
        rid = '_CCSM4_{}'.format(rcp)
        title_str = rcp[:-1] + '.' + rcp[4:] + ' - ' + '2016'
#        graphics.plot_modeloutput_map(gdir, title=title_str, filesuffix=rid, modelyr=2016, ax=axs[rn,0], vmax=350)
        title_str = rcp[:-1] + '.' + rcp[4:] + ' - ' + '2050'
#        graphics.plot_modeloutput_map(gdir, title=title_str, filesuffix=rid, modelyr=2050, ax=axs[rn,1], vmax=350)
        title_str = rcp[:-1] + '.' + rcp[4:] + ' - ' + '2100'
#        graphics.plot_modeloutput_map(gdir, title=title_str, filesuffix=rid, modelyr=2100, ax=axs[rn,2], vmax=350)
        rn += 1
#    plt.tight_layout()
#    plt.savefig(os.path.join(figfp, gdir.rgi_id + '_GCM_model_output_plot.png'), dpi=300)


# group by year and pass a sum or mean function to columns
rcp26_mean = rcp26_result.groupby('calendar_year').agg(
                                                    volume=('volume', 'sum'),
                                                    area=('area', 'sum'),
                                                    length=('length', 'mean')
                                                    )
rcp45_mean = rcp45_result.groupby('calendar_year').agg(
                                                    volume=('volume', 'sum'),
                                                    area=('area', 'sum'),
                                                    length=('length', 'mean')
                                                    )
rcp60_mean = rcp60_result.groupby('calendar_year').agg(
                                                    volume=('volume', 'sum'),
                                                    area=('area', 'sum'),
                                                    length=('length', 'mean')
                                                    )
rcp85_mean = rcp85_result.groupby('calendar_year').agg(
                                                    volume=('volume', 'sum'),
                                                    area=('area', 'sum'),
                                                    length=('length', 'mean')
                                                    )


def filt(df, col):
    """
    Pass rolling filter to dataframe column

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    col : str
        DESCRIPTION.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    df[col] = df[col].rolling(5, center=True, min_periods=1).mean()
    return df

# filter area and length columns, info https://oggm.org/tutorials/notebooks/area_length_filter.html
rcp26_mean = filt(rcp26_mean, 'area')
rcp26_mean = filt(rcp26_mean, 'length')
rcp45_mean = filt(rcp45_mean, 'area')
rcp45_mean = filt(rcp45_mean, 'length')
rcp60_mean = filt(rcp60_mean, 'area')
rcp60_mean = filt(rcp60_mean, 'length')
rcp85_mean = filt(rcp85_mean, 'area')
rcp85_mean = filt(rcp85_mean, 'length')

# merge dataframes for saving
merged = rcp26_mean.merge(rcp45_mean, how='outer', on='calendar_year', suffixes=(None, '_rcp45'))
merged = merged.merge(rcp60_mean, how='outer', on='calendar_year', suffixes=(None, '_rcp60'))
merged = merged.merge(rcp85_mean, how='outer', on='calendar_year', suffixes=(None, '_rcp85'))

# add suffix '_rcp26' to first columns
merged = merged.rename(columns={'volume': 'volume_rcp26', 'area': 'area_rcp26', 'length': 'length_rcp26'})


# save dataframe to file in case of later need
outname = os.path.join(os.path.dirname(figfp), 'rcp_result_sum.csv')
merged.to_csv(outname, sep=';')

# read csv file
rcp_result = pd.read_csv(outname, sep=';')

# convert cubic and square meters to cubic and square kilometers
for i in rcp_result.columns:
    if 'volume' in i:
        rcp_result[i] = rcp_result[i] / 10**9
    if 'area' in i:
        rcp_result[i] = rcp_result[i] / 10**6
        
# plot
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2,2)

ax1 = fig.add_subplot(gs[0,:])
ax1.plot(rcp_result.calendar_year, rcp_result.volume_rcp26, label='rcp2.6', color='b')
ax1.plot(rcp_result.calendar_year, rcp_result.volume_rcp45, label='rcp4.5', color='#006633')
ax1.plot(rcp_result.calendar_year, rcp_result.volume_rcp60, label='rcp6.0', color='#994C00')
ax1.plot(rcp_result.calendar_year, rcp_result.volume_rcp85, label='rcp8.5', color='k')
ax1.set_ylabel('volume (km$^3$)')
ax1.legend()

ax2 = fig.add_subplot(gs[1,0])
ax2.plot(rcp_result.calendar_year, rcp_result.area_rcp26, label='rcp2.6', color='b')
ax2.plot(rcp_result.calendar_year, rcp_result.area_rcp45, label='rcp4.5', color='#006633')
ax2.plot(rcp_result.calendar_year, rcp_result.area_rcp60, label='rcp6.0', color='#994C00')
ax2.plot(rcp_result.calendar_year, rcp_result.area_rcp85, label='rcp8.5', color='k')
ax2.set_ylabel('area (km$^2$)')

ax3 = fig.add_subplot(gs[1,1])
ax3.plot(rcp_result.calendar_year, rcp_result.length_rcp26, label='rcp2.6', color='b')
ax3.plot(rcp_result.calendar_year, rcp_result.length_rcp45, label='rcp4.5', color='#006633')
ax3.plot(rcp_result.calendar_year, rcp_result.length_rcp60, label='rcp6.0', color='#994C00')
ax3.plot(rcp_result.calendar_year, rcp_result.length_rcp85, label='rcp8.5', color='k')
ax3.set_ylabel('length (m)')

fig.suptitle('Projected changes 2016-2100')
plt.savefig(os.path.join(os.path.dirname(figfp), 'projected_change_km_2016_2100.png'), dpi=300)


# print some results

# set calendar year as index
rcp_result = rcp_result.set_index('calendar_year')

def printRowRatio(df, col1, col2, year1, year2):
    """
    Print percentage difference between selected years in selected columns

    Parameters
    ----------
    df : DataFrame
        DESCRIPTION.
    col1 : String
        column where to take start value.
    col2 : String
        column where to take end value.
    year1 : Integer
        Starting year to select
    year2 : Integer
        End year to select
        
    Returns
    -------
    Prints string of the result

    """

    start_vol = df.loc[[year1], [col1]]
    start_vol = start_vol[col1].iloc[0]
    end_vol = df.loc[[year2], [col2]]
    end_vol = end_vol[col2].iloc[0]
    
    percentage = round(100- (end_vol / start_vol * 100), 2)
    
    print('Start value: ' + str(round(start_vol, 2)))
    print('End value: ' + str(round(end_vol, 2)))
    print('Estimated change: ' + str(percentage) + ' %' ) 
    
    return percentage

printRowRatio(rcp_result, 'length_rcp26', 'length_rcp26', 2015, 2100)

#### get change between first and last row in percentage to table

# empty dataframe
perc_result = pd.DataFrame(columns=['volume', 'area', 'length'], index=['rcp26', 'rcp45', 'rcp60', 'rcp85'])

for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
        
    for i in rcp_result.columns:
        # if rcp found in column name, compute percentage and assign to df
        if rcp in i:
            # get result type from column name
            rtype = i[:-6]
            
            # subset and compute percentage between first and last row
            selcol = rcp_result[i]
            sel_percent = round(100- (selcol.iloc[-1] / selcol.iloc[0] * 100), 2)
    
            # assign result to certain row and column
            perc_result[rtype].loc[rcp] = sel_percent
        # return to the beginning if rcp not in column name
        else:
            continue

# save dataframe to file
outname2 = os.path.join(figfp, 'rcp_result_diff_percentage.csv')
perc_result.to_csv(outname2, sep=';')























