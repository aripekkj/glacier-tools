#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 13:45:16 2021

Use OGGM to compute centerlines and bedshape, then extract elevation 
from multiple DEMs to centerlines and plot

Notes:
    For a larger number of glaciers the script is a bit slow since there's quite many 
    modelling steps

@author: apj
"""


import os
import glob
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import oggm
import rasterio
import shapely
from shapely.geometry import LineString
import fiona
from shapely import geometry
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import gis
from oggm.core import gcm_climate, flowline
from pyproj import CRS
from functools import partial, wraps
from shapely.ops import transform as shp_trafo
from salem import wgs84
from rgitools.funcs import compute_intersects
import shutil

# fig fp
figfp = r'/Users/apj/Documents/_HY/Greenland/OGGM/figs'

# DEM folder
dem_dir = r'/Users/apj/Documents/_HY/Greenland/OGGM/*.tif'

cfg.initialize(logging_level='WORKFLOW')
cfg.PARAMS['border'] = 20

# read dem
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/study_area_dems_coregistered/80s_dem_studyarea_wgs84.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/coreg_vgridshift/1953_masked_2212_2016dem_wgs84_ellipsoid_utm_25m_glaciermasked_nuth_x+1.82_y-13.07_z-2.21_align.tif'
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/1953_dem_wgs84.tif'
custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/OGGM/1953_dem_interpolated_wgs84.tif'

custom_dem_path

# set custom dem path to PATHS, so it can be used in tasks. tasks.define_glacier_region parameter source='user' will read dem from cfg.PATHS['dem_file']
cfg.PATHS['dem_file'] = custom_dem_path

# do not use intersects
cfg.PARAMS['use_intersects'] = True
cfg.PARAMS['use_rgi_area'] = False

# read shapefile
#fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/1953_extent_utm_edited_final_edit_wgs.shp'
fp53 = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/1953_extent_utm_edited_final_edit_wgs.shp'

# filepath
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/final/*_final_edit_wgs*.shp'

# Create a dictionary of outlines
# new dictionary
shapes = dict()
# read files to dictionary
for outline in glob.glob(fp):
    key = os.path.basename(outline)[:4]
    shapes[key] = gpd.read_file(outline)
      
# function to get lowest point inside a masked polygon
def getLowestPoint(fp_raster, gdf, rgi_id):
    """
    Get lowest value of raster inside polygon mask

    Parameters
    ----------
    fp_raster : str
        Raster filepath.
    gdf : GeoDataFrame
        GeoDataFrame of polygons.
    rgi_id : str
        ID to select polygon.

    Returns
    -------
    masked_min : float
        minimum value of raster.

    """
    # select glacier by rgi_id
    sel = gdf[gdf.RGIId == rgi_id]

    with rasterio.open(fp_raster) as src:
        glac_mask, glac_out_transform = rasterio.mask.mask(src, sel.geometry, crop=False)
        glac_nodata = src.nodata

    masked = glac_mask[0]
    masked[(masked == glac_nodata)] = np.nan    
    
    # get lowest value
    masked_min = np.nanmin(masked)
    return masked_min
 
gdf = gpd.read_file(fp53) # read one of the outline files to use as 'default'

# exclude glaciers smaller than 0.1 km2
#gdf = entity[entity['Area'] >= 0.1]

# select or deselect ids
excludeIDs = ['RGI60-05.01920', 'RGI60-05.01987_1', 'RGI60-05.02213', 'RGI60-05.02328_1', 'RGI60-05.02280_1']
gdf = gdf[gdf.RGIId.isin(excludeIDs)]

#outline_year = os.path.basename(fp)[:4]

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('User_surf_change', reset=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

# empty GeoDataFrame to save flowlines
flowlines = gpd.GeoDataFrame()

for rgiid in gdf['RGIId']:
  
    # select glacier
    entity = gdf[gdf['RGIId'] == rgiid]
    
    # glacier directory
    gdirs = workflow.init_glacier_directories(entity)#, prepro_base_url=base_url) # <- from prepro_level uses preprocessed files
    
    gdir = gdirs[0]
    
    # glacier region
    tasks.define_glacier_region(gdir, source='USER') # set user dem
    
    # plot glacier
    #graphics.plot_domain(gdirs) 
    
    # print gdir files
    #print(os.listdir(gdirs.dir))
    
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
        try: # try statement allows to skip errors
            # The order matters!
            workflow.execute_entity_task(task, gdirs)
        except: # if exception is raised, add ID to list and return to beginning of loop
            excludeIDs.append(rgiid)
            pass
            continue
    #graphics.plot_centerlines(gdir, figsize=(8, 7), use_flowlines=True, add_downstream=True)
    
    ### Climate tasks ###
    # get tstars data to working dir
    tstar_url = 'https://cluster.klima.uni-bremen.de/~oggm/ref_mb_params/oggm_v1.4/RGIV62/CRU/centerlines/qc3/pcp2.5'
    workflow.download_ref_tstars(base_url=tstar_url)
    
    # run climate related entity tasks
    try: # try statement allows to skip errors
        workflow.climate_tasks(gdirs) # Downloads some files on the first time!
    except: # if exception is raised, add ID to list and return to beginning of loop
        excludeIDs.append(rgiid)
        pass
        continue
    
    ### Mass balance ###
    from oggm.core.massbalance import MultipleFlowlineMassBalance
    mbmod = MultipleFlowlineMassBalance(gdir, use_inversion_flowlines=True)
    
    years = np.arange(1953, 2016)
    mb_ts = mbmod.get_specific_mb(year=years)
    #plt.plot(years, mb_ts); plt.ylabel('SMB (mm yr$^{-1}$)')
    
    ### Ice thickness ###
    list_talks = [
             tasks.prepare_for_inversion,  # This is a preprocessing task
             tasks.mass_conservation_inversion,  # This does the actual job
             tasks.filter_inversion_output  # This smoothes the thicknesses at the tongue a little
             ]
    for task in list_talks:
        workflow.execute_entity_task(task, gdirs)
    
    # plot
    #graphics.plot_inversion(gdirs, figsize=(8, 7))
    
    
    # from tutorial
    tasks.init_present_time_glacier(gdir)
    tasks.run_constant_climate(gdir, nyears=100, y0=2000);
    
    fmod = flowline.FileModel(gdir.get_filepath('model_run'))
    fmod.run_until(0)
    #graphics.plot_modeloutput_map(gdir, model=fmod) # plot
    
    
    # get glacier flowline and bed elevation
    # Main flowline from FlowlineMassBalance
    fl = fmod.fls[-1]
    i, j = fl.line.xy  # xy flowline on grid
    lons, lats = gdir.grid.ij_to_crs(i, j, crs='EPSG:4326')  # to WGS84
    
    df_coords = pd.DataFrame(index=fl.dis_on_line*gdir.grid.dx)
    df_coords.index.name = 'Distance along flowline'
    df_coords['lon'] = lons
    df_coords['lat'] = lats
    df_coords['bed_elevation'] = fl.bed_h
        
    # convert df_coords to geodataframe
    gdf_coords = gpd.GeoDataFrame(df_coords, geometry=gpd.points_from_xy(df_coords.lon, df_coords.lat))
    # reset index
    gdf_coords = gdf_coords.reset_index()
    # Add RGI ID to column so gpd can later be grouped by ID
    gdf_coords['RGIId'] = gdir.rgi_id    
    
    # append gdf to flowlines gdf
    flowlines = flowlines.append(gdf_coords)
    
    # Now store a time varying array of ice thickness, surface elevation along this line
    df_thick = pd.DataFrame(index=df_coords.index)
    df_surf_h = pd.DataFrame(index=df_coords.index)
    df_bed_h = pd.DataFrame()
    for year in range(0, 101):
        fmod.run_until(year)
        fl = fmod.fls[-1]
        df_thick[year] = fl.thick
        df_surf_h[year] = fl.surface_h
    
    
    # plot simulated elevation change with bed elevation, from tutorial
    #f, ax = plt.subplots()
    #df_surf_h[[0, 50, 100]].plot(ax=ax);
    #df_coords['bed_elevation'].plot(ax=ax, color='k');
    #plt.title('Ice thickness');
    
    # loop through DEMs and extract flow line elevations
    for dem_fp in glob.glob(dem_dir):
        #with rasterio.open(dem_fp) as src:
        #    dem = src.read(1)
        dem = rasterio.open(dem_fp)
        #print('Processing DEM: ' + dem_fp)
        dem_colname = 'dem_' + os.path.basename(dem_fp)[:4]
        # extract DEM elevations
        for idx,rows in df_coords.iterrows():
            x = rows['lon'] 
            y = rows['lat']
            #print(x, y)
            
            # get row and column corresponding x y coordinate
            row, col = dem.index(x,y)    
            # extract DEM value
            dem_val = dem.read(1)[row,col]
            # assign dem value to column
            df_coords.loc[df_coords.index == idx, dem_colname] = dem_val

    # plot simulated elevation change with bed elevation
    f, ax = plt.subplots()
    #df_surf_h[0].plot(ax=ax, label='Simulated year 0', ls='--', linewidth= 0.9, color='#B2BABB')
    #df_surf_h[32].plot(ax=ax, label='Simulated year 32', ls='--', linewidth= 0.9, color='#707B7C')
    #df_surf_h[63].plot(ax=ax, label='Simulated year 63', ls='--', linewidth= 0.9, color='#626567')
    df_coords['dem_1953'].plot(ax=ax, linewidth= 1, color='k')
    df_coords['dem_1985'].plot(ax=ax, linewidth= 1, color='#994C00')
    df_coords['dem_2016'].plot(ax=ax, linewidth= 1, color='#006633')
    df_coords['bed_elevation'].plot(ax=ax, ls='--', color='#707B7C')
    plt.ylabel('Altitude (m)')
    plt.legend()
    plt.title('Ice thickness change on ' + gdir.rgi_id)
    #plt.savefig(os.path.join(figfp, '_' + gdir.rgi_id + '_l_icethickness.png'), dpi=150)    

# group flowline gdf by RGIId and convert Point geometry to LineString
flowlines_grouped = flowlines.groupby(['RGIId'])['geometry'].apply(lambda x: LineString(x.tolist()))
gdf_flowlines = gpd.GeoDataFrame(flowlines_grouped, geometry='geometry')

# write GeoDataFrame to shapefile
outname = os.path.join(os.path.dirname(figfp), 'flowlines.shp')
gdf_flowlines.to_file(outname, driver='ESRI Shapefile')













