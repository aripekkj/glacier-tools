# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:15:51 2020

Create a folium map with popup images


Updates:
    Added MarkerCluster class with popups 15.6.2020

To Do:
    Add centerlines to map

@author: Ap
"""

import base64
import folium
import pandas as pd
import geopandas as gpd
from folium import IFrame
from folium.plugins import MarkerCluster
import glob
import os
import numpy as np

# function to extract row by string
def getPointByString(String):
    row = f.loc[f['RGIId'] == String]
    row_point = row.centroid
    return [row, row_point];
    
# polygon filepath
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/05_rgi60_nuussuaq_wgs84_studyarea.shp'
# centerlines filepath
#line_fp = r'/Users/apj/Documents/_HY/Greenland/'

# read file
f = gpd.read_file(fp)

# select only columns that are needed
#fire = firecount[['WDPA_PID', 'NAME', 'diff_2020_', 'cat', 'geometry']]
#fire['geo_id'] = fire.index.astype(str)
# create columns for x and y coordinates
f['x'] = f.geometry.centroid.x
f['y'] = f.geometry.centroid.y

# map location
loc = 70.00, -52.50

# set path for glob to browse through files
img_path = r'/Users/apj/Documents/_HY/Greenland/centerlines/figures/*.png' # look only files that end with _s.png

# set categories for choropleth map legend
#bins = list([0, 25, 50, max(fire['diff20-19'])]) # not working

# folium map
m = folium.Map(location=loc, zoom_start=7, tiles='Stamen Terrain')


#add polygons to map
folium.Choropleth(
    geo_data=f,
    data=f,
    columns=['RGIId', 'Surging'],
    key_on='feature.id',
    fill_color='YlOrRd',
    fill_opacity=0.7,
    line_opacity=0.5,
    legend_name='RGI surging status',
    bins=4,
    reset=True
).add_to(m)


# Convert polygon to GeoJson and add as a transparent layer to the map with tooltip
folium.features.GeoJson(f,
                        name='Labels',
                        style_function=lambda x: {'color':'transparent','fillColor':'transparent','weight':0},
                        tooltip=folium.features.GeoJsonTooltip(fields=['RGIId', 'Surging'],
                                                                aliases = ['Glacier ID', 'Surging status'],
                                                                labels=True,
                                                                sticky=False
                                                                            )
                       ).add_to(m)




# empty geodataframe for point coordinate rows
point_coords = gpd.GeoDataFrame()
# empty lists to store objects
popuplist = []
iconlist = []
# browse through files in folder
for filename in glob.glob(img_path): 
    file_only = os.path.basename(filename[:-4]) # Get filename only to extract row
    img_fp = os.path.abspath(filename) # get full filepath to link image
    
    result = getPointByString(file_only) # use function to extract row and get point coordinates
    marker_row = result[0]
    marker_point = result[1]
    
    # add row to dataframe
    point_coords = point_coords.append(marker_row)
    
    # open png image
    encoded = base64.b64encode(open(img_fp, 'rb').read())

    # add marker and png image as popup
    html = '<img src="data:image/png;base64,{}">'.format
    iframe = IFrame(html(encoded.decode('UTF-8')), width=620, height=420)
    popup = folium.Popup(iframe, max_width= 600)
    popuplist.append(popup)
    icon = folium.Icon(color="red", icon="info-sign")
    iconlist.append(icon)
    # simple marker
    #marker = folium.Marker(location=[marker_point.y, marker_point.x],
    #                       popup=popup,
    #                       icon=icon) # 
    #marker.add_to(m)
    

# create a list of coordinates from extracted rows
locations = list(zip(point_coords["y"], point_coords["x"]))

# create marker_cluster class and add it to map
marker_cluster = MarkerCluster(locations, popups=popuplist, icons=iconlist)
m.add_child(marker_cluster)

# save map
m.save("folium_map_glacier.html")








