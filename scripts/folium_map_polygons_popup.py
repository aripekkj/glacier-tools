# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:15:51 2020

Create a folium map with popup images

Notes: 
    Use lat, lon coordinates in popup markers.

Updates:
    Script can now add lines to map 23.9.2020
    Added MarkerCluster class with popups 15.6.2020
    

To Do:
    Set popup marker locations somewhere on LineString instead of line centroids

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
def getPointByString(DataFrame, String):
    row = DataFrame.loc[DataFrame['RGIID2'] == String]
    row_point = row.centroid
    return [row, row_point];
  
# polygon filepath
fp = r'/Users/apj/Documents/_HY/Greenland/outlines/edited_glacier_divides/'
# centerlines filepath
line_fp = r'/Users/apj/Documents/_HY/Greenland/centerlines/edited_glacier_divides/unique_ids/cl_50s_wgs84_edit_newUniqueIDs.shp'

# read files
poly = gpd.read_file(fp)
lines = gpd.read_file(line_fp)

# new column to create unique ID's also for tributary glaciers
#lines['RGIID_OFID'] = lines['RGIID'] + '_' + lines['ORIG_FID'].astype(str)

# create columns for x and y coordinates
lines['x'] = lines.geometry.centroid.x
lines['y'] = lines.geometry.centroid.y


# map location
loc = 70.50, -52.50

# set path for glob to browse through files
img_path = r'/Users/apj/Documents/_HY/Greenland/centerlines/figures/edited_divides_w_tributaries/*.png' # look only files that end with _s.png

# set categories for choropleth map legend
#bins = list([0, 25, 50, max(fire['diff20-19'])]) # not working

# folium map
m = folium.Map(location=loc, zoom_start=8, tiles='Stamen Terrain')


#add polygons to map
folium.Choropleth(
    geo_data=poly,
    data=poly,
    columns=['RGIId', 'Surging'],
    key_on='feature.properties.RGIId',
    fill_color='YlOrRd',
    fill_opacity=0.7,
    line_opacity=0.5,
    legend_name='RGI surging status',
    bins=4,
    reset=True
).add_to(m)


# Convert polygon to GeoJson and add as a transparent layer to the map with tooltip
folium.features.GeoJson(poly,
                        name='Labels',
                        style_function=lambda x: {'color':'transparent','fillColor':'YlOrRd','weight':1},
                        tooltip=folium.features.GeoJsonTooltip(fields=['Surging'],
                                                                aliases = ['Surging status'],
                                                                labels=True,
                                                                sticky=False
                                                                            )
                       ).add_to(m)


# Convert lines to GeoJson and add as a transparent layer to the map with tooltip
folium.features.GeoJson(lines,
                        name='Lines',
                        style_function=lambda x: {'color':'blue','weight':2},
                        tooltip=folium.features.GeoJsonTooltip(fields=['RGIID'],
                                                                aliases = ['Glacier ID'],
                                                                labels=True,
                                                                sticky=False
                                                                            )
                       ).add_to(m)



# This part is to add credit text on the map
"""
html = '<div style="position: fixed; bottom: 30px; right: 5px; width: 200px; height: 60px; \
    background-color: #FFFFFF00; z-index:9000; line-height: 10px"> \
        <font size="1">\
        Data: \
        <br>VIIRS Active Fire product \
        <br><a href="https://earthdata.nasa.gov/firms">https://earthdata.nasa.gov/firms</a> \
        <a href="https://earthdata.nasa.gov/earth-observation-data/near-real-time/firms/v1-vnp14imgt">DOI</a>\
        <br>UNEP-WCMC and IUCN (2020) \
        <br><a href="https://www.protectedplanet.net">www.protectedplanet.net</a> \
        <br>Data visualization: Ari-Pekka Jokinen\
        </font>\
        </div>' 
m.get_root().html.add_child(folium.Element(html))
"""


# empty geodataframe for point coordinate rows
point_coords = gpd.GeoDataFrame()
# empty lists to store objects
popuplist = []
iconlist = []
# browse through files in folder
for filename in glob.glob(img_path): 
    file_only = os.path.basename(filename[:-4]) # Get filename only to extract row
    img_fp = os.path.abspath(filename) # get full filepath to link image
    
    result = getPointByString(lines,file_only) # use function to extract row and get point coordinates
    marker_row = result[0]
    marker_point = result[1]
    
    # add row to dataframe
    point_coords = point_coords.append(marker_row)
    
    # open png image
    encoded = base64.b64encode(open(img_fp, 'rb').read())

    # add marker and png image as popup
    html = '<img src="data:image/png;base64,{}">'.format
    iframe = IFrame(html(encoded.decode('UTF-8')), width=450, height=300)
    popup = folium.Popup(iframe, max_width= 450)
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
m.save("folium_map_glacierlines_3010.html")








