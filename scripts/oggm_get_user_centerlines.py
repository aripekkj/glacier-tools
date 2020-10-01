
"""
Created on Tue Sep 29 17:47:44 2020


Notes:
    1.10. 2020
    Creates centerlines on user imported shape. Not sure about DEM
    
    1.10.2020
    Seems to work with python 3.6 and latest OGGM build
    
    29.9.2020
    When usign the latest build of OGGM and python 3.7:
    
    "PicklingError: Can't pickle <class 'oggm.utils._workflow.GlacierDirectory'>: it's not the same object as oggm.utils._workflow.GlacierDirectory"
    
    

@author: apj
"""

import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import oggm
from oggm import cfg, utils, workflow, tasks, graphics
from oggm.core import gis
from pyproj import CRS

from functools import partial, wraps
from shapely.ops import transform as shp_trafo
from salem import wgs84


cfg.initialize(logging_level='WORKFLOW')
cfg.PARAMS['border'] = 10


# out fp
out_fp = r'/Users/apj/Documents/_HY/Greenland/centerlines/centerlines_50s_edited.shp'

# read dem
custom_dem_path = r'/Users/apj/Documents/_HY/Greenland/demcoreg_result/study_area_dems_coregistered/50s_dem_wgs84.tif'
custom_dem_path

cfg.PATHS['dem_file'] = custom_dem_path

# do not use intersects
cfg.PARAMS['use_intersects'] = False

# read shapefile
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/05_rgi60_studyarea_50s_extent_wgs84_edited.shp'
entity = gpd.read_file(fp)

# test with one glacier
#entity = entity[entity['RGIId'] == 'RGI60-05.02208_2']

cfg.PATHS['working_dir'] = utils.gettempdir('user', reset=True)
gdirs = workflow.init_glacier_directories(entity)


# tasks to be executed
list_tasks = [
         tasks.define_glacier_region,
         tasks.glacier_masks, 
         tasks.compute_centerlines,
         tasks.initialize_flowlines,
         tasks.catchment_area,
         tasks.catchment_width_geom,
         tasks.catchment_width_correction,
         tasks.compute_downstream_line,
         ]
for task in list_tasks:
    # The order matters!
    workflow.execute_entity_task(task, gdirs)


# plot. suitable for plotting one glacier
#f, ax = plt.subplots()
#da_user = xr.open_rasterio(gdirs.get_filepath('dem'))
#da_user.plot(cmap='terrain', ax=ax);
#gdirs.read_shapefile('outlines').plot(ax=ax, color='none', edgecolor='black');
#graphics.plot_centerlines(gdirs, use_flowlines=False)



# get centerlines
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















