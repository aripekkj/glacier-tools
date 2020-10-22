#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:23:56 2020

Create numbering for tributary glaciers


@author: apj
"""

# input file
fp_line = r'/Users/apj/Documents/_HY/Greenland/centerlines/edited_glacier_divides/centerlines_50s_edited_utm_11102020.gpkg'

# output file
lineout = r'/Users/apj/Documents/_HY/Greenland/centerlines/edited_glacier_divides/centerlines_50s_utm_tributaryIDs_22102020.gpkg'

# read lines
lines = gpd.read_file(fp_line)

# change 'MAIN' column to object type
lines = lines.astype({"MAIN": object})
lines.dtypes

# set running number
running_number = 1
# If glacier has more than one flowline, number lines in 'MAIN' column 1.1, 1.2, 1.3, etc.
for row in lines.itertuples():
    # count how many lines are found with the RGIID
    linecount = lines.RGIID.value_counts()[row.RGIID]
    print('RGIID: ' + str(row.RGIID) + ', Line count: ' +  str(linecount))
    if linecount > 1:
        lines.at[row.Index, 'MAIN'] = str(1) + '.' + str(running_number) 
        print(str(1) + '.' + str(running_number))    
        if running_number == linecount:
            running_number == 1
        else:
            running_number += 1
        continue
    
    running_number = 1
    
# write lines to file
lines.to_file(lineout, driver='GPKG')
