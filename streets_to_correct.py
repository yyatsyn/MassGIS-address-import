# This file prepare a list of streets which names do not match with
# .. names in MassGIS; Later those streets can be used for a MapRoulette challenge
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import re
from shapely.geometry import Point, LineString, Polygon
import difflib

# Some list and functions used later for fuzzy matches 
to_remove = ["-", ",", "\.", "'", " "]
# most common words from str_all using from collections import Counter; Counter(test.split()).most_common()[0:200]
streets = ["Street", "Road", "Avenue", "Drive", "Lane", "Way", "Circle",
            "Trail", "Court", "Place", "Highway", "Path", "Turnpike", 
            "Boulevard", "Parkway", "Drive", "Driveway", "Extension", 
            "Loop", "Square", "Alley", "Terrace"] 
streets = [" " + x.upper() + "$" for x in streets]
streets.extend(to_remove)
to_remove = streets
# most common abbreviations
abbrevs = { " HWY$": " HIGHWAY", " HWY.$": " HIGHWAY",
            " ST$": " STREET", " ST.$": " STREET",
            " STR$": " STREET", " STR.$": " STREET",
            " RD$": " ROAD", " RD.$": " ROAD",
            " AVE$": " AVENUE", " AVE.$": " AVENUE",
            " LN$": " LANE", " LN.$": " LANE",
            " DR$": " DRIVE", " DR.$": " DRIVE",
            " TR$": " TRAIL", " TR.$": " TRAIL",
            " BLVD$": " BOULEVARD", " BLVD.$": " BOULEVARD",
            " EXT$": " EXTENSION", " EXT.$": " EXTENSION",
            " SQ$": " SQUARE", " SQ.$": " SQUARE",
            " PL$": " PLACE", " PL.$": " PLACE"}

def remove_simple_mismatches(street_name, to_remove, abbrevs):
    if street_name != None:
        street_name = street_name.upper()
        for k, v in abbrevs.items():
            street_name = re.sub(k, v, street_name)
        
        for i in to_remove: 
            street_name = re.sub(i, "", street_name)
            #print(i)
            #print(street_name)
    else:
        street_name = ""
    return street_name

#-----------------------------------------------------------------------
str_all = gpd.read_file(path+'all_streets_ma.shp')
cur = gpd.read_file(path + "cty_towns/problem_mgis_case_insensitive/all_combined.shp")

# number of instances with the same addr:street
cur["addr_n"] = cur.groupby(['cty_town','UL_STREET_'])['MASTER_ADD'].transform(len)
cur["x"] = cur.apply(lambda cur: cur["geometry"].xy[0][0], axis = 1)
cur["y"] = cur.apply(lambda cur: cur["geometry"].xy[1][0], axis = 1)
cur['point_id'] = cur.groupby(['x','y']).ngroup()
cur = cur.drop(['x','y'], axis=1)
cur["pnts_n"] = cur.groupby(['cty_town','UL_STREET_'])['point_id'].transform(lambda x: x.nunique())

drop_from_import = cur.loc[cur['addr_n']==1] # points with 1 address of a kind is usually erroneous
cur['geometry'] = cur['geometry'].apply(lambda x: x.coords[0])
cur_pnt = cur.loc[(cur['pnts_n']==1) & (cur["addr_n"]>1)]
cur_lns = cur.loc[cur['pnts_n']==2]
cur_plg = cur.loc[cur['pnts_n']>2]

to_check_plg = gpd.GeoDataFrame(cur_plg.groupby(['cty_town','UL_STREET_'])['geometry']\
    .apply(lambda x: Polygon(x.tolist()))).reset_index()
to_check_plg['geometry'] = to_check_plg['geometry'].convex_hull
to_check_lns = gpd.GeoDataFrame(cur_lns.groupby(['cty_town','UL_STREET_'])['geometry']\
    .apply(lambda x: LineString(x.tolist()))).reset_index()
to_check_pnt = gpd.GeoDataFrame(cur_pnt.groupby(['cty_town','UL_STREET_'])['geometry']\
    .apply(lambda x: Point(x.tolist()))).reset_index()

to_check = to_check_lns.append(to_check_pnt, ignore_index=True)
to_check['geometry'] = to_check['geometry'].buffer(0.001)
to_check_plg['geometry'] = to_check_plg['geometry'].buffer(0.001)
to_check = to_check.append(to_check_plg, ignore_index=True)
to_check['long'] = to_check['geometry'].length

# sjoin OSM streets and MassGIS-based polygons with fuzzy-unmatched addresses
str_to_check = gpd.sjoin(str_all, to_check, how="inner", op="intersects")

# match MassGIS and OSM streets by names
str_to_check['name'] = np.where(str_to_check['name'].isna(), "", str_to_check['name'])
str_to_check['name_osm'] = str_to_check["name"]\
    .apply(lambda x: remove_simple_mismatches(x, to_remove, abbrevs))
str_to_check['name_mgis'] = str_to_check["UL_STREET_"]\
    .apply(lambda x: remove_simple_mismatches(x, to_remove, abbrevs))
str_to_check['sim_name'] = str_to_check\
    .apply(lambda x: difflib.SequenceMatcher(None, x['name_osm'], x['name_mgis']).ratio(), axis=1)

# keep streets with closest match in names
str_to_check["sim_max"] = str_to_check.groupby(['cty_town','UL_STREET_'])['sim_name'].transform(max)
str_to_check2 = str_to_check.loc[str_to_check["sim_max"]==str_to_check["sim_name"]]
str_to_check3 = str_to_check.loc[(str_to_check["sim_max"]==str_to_check["sim_name"]) & (str_to_check["sim_max"]>0.8)]
#xxx=str_to_check2.loc[str_to_check2['cty_town']=='MIDDLESEX, BILLERICA', ['name', 'UL_STREET_', 'full_id', 'cty_town']]
#str_to_check2[['name', 'UL_STREET_', 'full_id', 'cty_town']]
#xxx
#list(xxx.iloc[:,2])
str_to_check2.drop(['geometry'], axis=1).to_csv(path + "mismatched_streets_case_insensitive_raw.csv", index=False, sep = ';')
str_to_check3.to_csv(path + "mismatched_streets_case_insensitive.csv", index=False)

# aggregated to find a street that intersects the convex hull with the longest segment
#str_to_check['long2'] = str_to_check['geometry'].length
#str_to_check['long2'] = str_to_check.groupby(['cty_town','UL_STREET_','name'])['long2'].transform(sum)
#str_to_check["len_max"] = str_to_check.groupby(['cty_town','UL_STREET_'])['long2'].transform(max)
#str_to_check = str_to_check.loc[str_to_check["len_max"]==str_to_check["long2"]]
tojson = str_all.loc[str_all['full_id'].isin(str_to_check3['full_id'])]

# import matplotlib.pyplot as plt
# df = to_check.loc[to_check['long']>0.1]
# df.plot()
# plt.show()
