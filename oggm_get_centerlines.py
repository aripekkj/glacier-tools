#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:01:49 2019

Extract OGGM centerlines. Script based on OGGM tutorial


@author: apj
"""

# import packages
import gdal
import oggm
import numpy as np
from oggm import cfg
import geopandas as gpd
from functools import partial, wraps
from shapely.ops import transform as shp_trafo
from salem import wgs84
from pyproj import CRS

# create configuration file
cfg.initialize()

# locate configuration file
cfg.CONFIG_FILE

# out fp
out_fp = r'/Users/apj/Documents/_HY/Greenland/centerlines/centerlines_rgi.shp'

# test run. Read the parameter file and make it available for tools.
from oggm import cfg, utils
cfg.initialize(logging_level='WORKFLOW')

# some parameters
cfg.PARAMS['continue_on_error'] = True
cfg.PARAMS['border'] = 10

# define user dem
#custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/AERODEM/aerodem_1985_wgs84.tif'
#custom_dem_path

#cfg.PATHS['dem_file'] = custom_dem_path


#cfg.PATHS['working_dir'] = utils.gettempdir('user')

""" Workflow """
from oggm import workflow

# set working directory
#cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-GettingStarted', reset=True)
cfg.PATHS['working_dir'] = '/Users/apj/OGGM/'
cfg.PATHS['working_dir']

# Define glaciers for the run. These are called by Randolph glacier inventory id
#rgi_ids = ['RGI60-11.01328', 'RGI60-11.00897']
#rgi_ids = ['RGI60-05.01125']




######### Select glaciers 
# 1 read shapefile
fp = r'/Users/apj/Documents/_HY/Greenland/masks/nuussuaq_rgi_wgs84.shp'
shape = gpd.read_file(fp)
shape.columns
shape['RGIId'].head()
#shape['DEM_SOURCE'] = 'USER'
# assign RGI Ids to list



# recompute glacier area. Important when using edited rgi!
cfg.PARAMS['use_rgi_area'] = False
# This is the default anyway, but we set it here to be sure
cfg.PARAMS['use_intersects'] = True

# Glacier directories from_prepro_level=1, prepro_border=80, prepro_rgi_version='61' from_tar=False
gdirs = workflow.init_glacier_regions(shape)

# look contents of gdirs
type(gdirs)

# access item on the gdirs list
gdir = gdirs[0]  # 
print('Path to the DEM:', gdir.get_filepath('dem'))

# look at the attributes of the glacier
gdir
gdir.rgi_date # date at which the outlines are valid

# plot the glacier location and outline
from oggm import graphics
#graphics.plot_googlemap(gdir, figsize=(8, 7))

""" Tasks """
from oggm import tasks

# run the glacier_masks task on all gdirs
workflow.execute_entity_task(tasks.glacier_masks, gdirs);

# the command wrote a new file in our glacier directory, providing raster masks of the glaciers
print('Path to the masks:', gdir.get_filepath('gridded_data'))

# It is also possible to apply several tasks sequentially (i.e. one after an other) on our glacier list:
list_talks = [
         tasks.compute_centerlines,
         tasks.initialize_flowlines,
         #tasks.compute_downstream_line,
         ]
for task in list_talks:
    # The order matters!
    workflow.execute_entity_task(task, gdirs)

# plot the results
#graphics.plot_centerlines(gdirs, figsize=(8, 7), use_flowlines=True, add_downstream=False)

# the glacier directories now have more files in them. Let's look
import os
print(os.listdir(gdir.dir))

olist = []
def _get_centerline_lonlat(gdir):
    
    # quick and dirty solution to export centerlines to shapefiles
    cls = gdir.read_pickle('inversion_flowlines')
    
    for j, cl in enumerate(cls[::-1]):
            mm = 1 if j == 0 else 0
            gs = gpd.GeoSeries()
            gs['RGIID'] = gdir.rgi_id
            gs['LE_SEGMENT'] = np.rint(np.max(cl.dis_on_line) * gdir.grid.dx)
            gs['MAIN'] = mm
            tra_func = partial(gdir.grid.ij_to_crs, crs=wgs84)
            gs['geometry'] = shp_trafo(tra_func, cl.line)
            olist.append(gs)
            
    return olist

for gdr in gdirs:
    _get_centerline_lonlat(gdr)


# assign olist to Geodataframe 
cl_gdf = gpd.GeoDataFrame(olist)

# get crs in well known text
def getWKT(epsg_code):
    crs = CRS.from_user_input(epsg_code)
    crs_wkt = crs.to_wkt(pretty=True)
    return crs_wkt

wgs84_wkt = getWKT(4326)

# set crs
cl_gdf.crs = wgs84_wkt

# write GeoDataFrame to Shapefile
cl_gdf.to_file(out_fp, driver='ESRI Shapefile')


# Another solution to write centerlines to shapefile (CRSError: Invalid input to create CRS: {'init': 'epsg:4326'})
#oggm.utils.write_centerlines_to_shape(gdirs, filesuffix='80s', path=True)











