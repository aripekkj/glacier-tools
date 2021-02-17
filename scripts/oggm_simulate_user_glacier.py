#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 07:43:00 2021

Run OGGM simulation for single glacier. Based on OGGM tutorials 'Getting started' and 'Run OGGM with GCM data'

Notes:
    13.2.2021
    OGGM v1.4 latest build works on python 3.7 (MacOs)
    
    13.2.2021
    workflow.download_ref_tstars(base_url=params_url) <- workflow.py required editing
    Added parameter follow_symlinks=True 
    
    13.2.2021
    Datasets in WGS84 (EPSG:4326)
    
To Do:
    compute glacier intersects for own inventory    
    Get climate data
    
@author: apj
"""

import os
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
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
figfp = r'/Users/apj/OGGM/figs/'

cfg.initialize(logging_level='WORKFLOW')
cfg.PARAMS['border'] = 10

# read dem
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/study_area_dems_coregistered/80s_dem_studyarea_wgs84.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1953_masked_2212_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+1.82_y-13.07_z-2.21_align.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/1953_dem_wgs84.tif'
custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/2016dem_wgs84_studyarea.tif'

custom_dem_path

# set custom dem path to PATHS, so it can be used in tasks. tasks.define_glacier_region parameter source='user' will read dem from cfg.PATHS['dem_file']
cfg.PATHS['dem_file'] = custom_dem_path

# do not use intersects
cfg.PARAMS['use_intersects'] = True
cfg.PARAMS['use_rgi_area'] = False

# read shapefile
#fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/1953_extent_utm_edited_final_edit_wgs.shp'
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/2016_extent_updated_divides_utm_21102020_fixed_final_edit_wgs84.shp'

entity = gpd.read_file(fp)

outline_year = os.path.basename(fp)[:4]

# check geometries
for idx, row in entity.iterrows():
    if entity.geometry.iloc[idx].type != 'Polygon':
        print(row['RGIId'] + row.geometry.type)

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('User_gcm_run', reset=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

#cfg.PATHS['working_dir'] = utils.gettempdir('user', reset=True)
#cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-GettingStarted', reset=True)

# test with one glacier
entity = entity[entity['RGIId'] == 'RGI60-05.02114']

# rgi_ids for default data
#rgi_ids = ['RGI60-05.01916']

# glacier directory
#base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/L3-L5_files/CRU/centerlines/qc3/pcp2.5/no_match/'
gdirs = workflow.init_glacier_directories(entity)#, prepro_base_url=base_url) # <- from prepro_level uses preprocessed files
#gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3, prepro_border=80) # this is loading default dataset
gdir = gdirs[0]
rgi_id = gdir.rgi_id
gdir.rgi_date = outline_year

# glacier region
tasks.define_glacier_region(gdirs[0], source='USER') # set user dem

# plot glacier
graphics.plot_domain(gdirs[0]) 

# print gdir files
print(os.listdir(gdirs[0].dir))

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
    # The order matters!
    workflow.execute_entity_task(task, gdirs)

# plot catchment areas
graphics.plot_catchment_areas(gdirs, figsize=(8, 7))


### Climate tasks ###
# get tstars data to working dir
tstar_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'
workflow.download_ref_tstars(base_url=tstar_url)

# run climate related entity tasks
workflow.climate_tasks(gdirs) # Downloads some files on the first time!

""" 

Tutorial stuff

### Mass balance calibration from tutorial ###

## This follows getting started tutorial and uses preprocessed data
#fpath = gdirs[0].get_filepath('climate_historical')
#ds = xr.open_dataset(fpath)
# Data is in hydrological years
# -> let's just ignore the first and last calendar years
#ds.temp.resample(time='AS').mean()[1:-1].plot();

# Fetch the reference t* list and associated model parameters
#params_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'
#workflow.download_ref_tstars(base_url=params_url)
# Now calibrate
#workflow.execute_entity_task(tasks.local_t_star, gdirs);
#workflow.execute_entity_task(tasks.mu_star_calibration, gdirs);

### Mass balance ###
from oggm.core.massbalance import MultipleFlowlineMassBalance
mbmod = MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)

years = np.arange(1953, 2016)
mb_ts = mbmod.get_specific_mb(year=years)
plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)')


### Ice thickness ###
list_talks = [
         tasks.prepare_for_inversion,  # This is a preprocessing task
         tasks.mass_conservation_inversion,  # This does the actual job
         tasks.filter_inversion_output  # This smoothes the thicknesses at the tongue a little
         ]
for task in list_talks:
    workflow.execute_entity_task(task, gdirs)

# plot
graphics.plot_inversion(gdir, figsize=(8, 7))

cfg.PARAMS['inversion_glen_a']

a_factor = np.linspace(0.1, 10., 100)
volume = []
for f in a_factor:
    # Recompute the volume without overwriting the previous computations
    v = tasks.mass_conservation_inversion(gdir, glen_a=f * cfg.PARAMS['inversion_glen_a'], write=False)
    volume.append(v * 1e-9)
plt.plot(a_factor, volume); plt.title('Total volume');
plt.ylabel('Volume (km$^3$)'); plt.xlabel('Glen A factor (1 = default)'); 


### Simulation ###
# Convert the flowlines to a "glacier" for the ice dynamics module
workflow.execute_entity_task(tasks.init_present_time_glacier, gdirs)

# run random simulation from default CRU data. Randomize climate for nyears
#workflow.execute_entity_task(tasks.run_random_climate, gdirs, nyears=200,
#                             y0=1953, output_filesuffix='_2000')


# GCM filepaths
bp = 'https://cluster.klima.uni-bremen.de/~oggm/cmip5-ng/pr/pr_mon_CCSM4_{}_r1i1p1_g025.nc'
bt = 'https://cluster.klima.uni-bremen.de/~oggm/cmip5-ng/tas/tas_mon_CCSM4_{}_r1i1p1_g025.nc'


# process gcm data
workflow.execute_entity_task(tasks.process_cmip_data, gdirs, fpath_precip=bp, fpath_temp=bt)

# run simulation from climate data
workflow.execute_entity_task(tasks.run_from_climate_data, gdirs, ys=1953, ye=2016,
                             output_filesuffix='_hist')


ds2000 = utils.compile_run_output(gdirs, filesuffix='_hist')

# plot volume and length modelling
f, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
ds2000.volume.plot.line(ax=ax1, hue='rgi_id');
ds2000.length.plot.line(ax=ax2, hue='rgi_id');

# model output plot
f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))
graphics.plot_modeloutput_map(gdir, filesuffix='_hist', modelyr=0, ax=ax1, vmax=350)
graphics.plot_modeloutput_map(gdir, filesuffix='_hist', modelyr=30, ax=ax2, vmax=350)
graphics.plot_modeloutput_map(gdir, filesuffix='_hist', modelyr=50, ax=ax3, vmax=350)
plt.tight_layout();
#plt.savefig(os.path.join(figfp, 'model_output' + rgi_id + '.png'), dpi=150)


### Sensitivity to temperature change ###

# repeat simulation with 0-5 and -0.5 temp change
workflow.execute_entity_task(tasks.run_random_climate, gdirs, nyears=200,
                             temperature_bias=0.5,
                             y0=1953, output_filesuffix='_p05');
workflow.execute_entity_task(tasks.run_random_climate, gdirs, nyears=200,
                             temperature_bias=-0.5,
                             y0=1953, output_filesuffix='_m05');

# compile run
dsp = utils.compile_run_output(gdirs, filesuffix='_p05')
dsm = utils.compile_run_output(gdirs, filesuffix='_m05')

# plot
f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))
ds2000.sel(rgi_id=rgi_id).volume.plot.line(ax=ax1, hue='rgi_id', label='Commitment');
ds2000.sel(rgi_id=rgi_id).area.plot.line(ax=ax2, hue='rgi_id');
ds2000.sel(rgi_id=rgi_id).length.plot.line(ax=ax3, hue='rgi_id');
dsp.sel(rgi_id=rgi_id).volume.plot.line(ax=ax1, hue='rgi_id', label='$+$ 0.5°C');
dsp.sel(rgi_id=rgi_id).area.plot.line(ax=ax2, hue='rgi_id');
dsp.sel(rgi_id=rgi_id).length.plot.line(ax=ax3, hue='rgi_id');
dsm.sel(rgi_id=rgi_id).volume.plot.line(ax=ax1, hue='rgi_id', label='$-$ 0.5°C');
dsm.sel(rgi_id=rgi_id).area.plot.line(ax=ax2, hue='rgi_id');
dsm.sel(rgi_id=rgi_id).length.plot.line(ax=ax3, hue='rgi_id');
ax1.legend();

"""


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
                                 #init_model_filesuffix='_historical',  # this is important! Start from 2020 glacier
                                 output_filesuffix=rid,  # recognize the run for later
                                );

# plot modelling results
f, ax1 = plt.subplots(1, 1, figsize=(14, 4))
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
    rid = '_CCSM4_{}'.format(rcp)
    ds = utils.compile_run_output(gdirs, input_filesuffix=rid)
    ds.isel(rgi_id=0).volume.plot(ax=ax1, label=rcp);
    #ds.isel(rgi_id=1).volume.plot(ax=ax2, label=rcp);
plt.legend();











