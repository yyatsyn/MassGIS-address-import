# Find MassGIS addresses that may contain errors
# (in addition to the massgis_osm_mismatches above)
# Use two approaches (as they cover most of erroneous MassGIS addresses that I've encountered)
# 1) 1-3 addresses with a given street name in a town/suburb,
# ... but no such street name in OSM (very likely to include correct MassGIS address, but missing OSM street)
# 2) address with a given street name in a town/suburb that is VERY far
# ... from other addresses with such street name in a given town/suburb,
# ... distance to the nearest 'neighbor' is 10 times larger that the median distance
# ... to the nearest neighbor among all addresses with a given street name in a given town/suburb
# This step returns shp-files with address points which need more thorough check
# After checking these addresses add the list of MASTER_ADD to the
# Google spreadsheet ("Street names mismatches: MassGIS vs OSM")

import numpy as np
import geopandas as gpd
import pandas as pd
import sys
import os
import re
from titlecase import titlecase
from set_path import path
from natsort import natsorted, ns # for sorting strings that have mixed numbers and letters
from itertools import compress, groupby, count
from shapely.geometry import Point
from fiona.crs import from_epsg

# read data on counties
ctys = pd.read_csv(path+"town_cty_id.csv")
cty_towns = list(ctys["cty_town"].unique())
ctys = list(ctys["cty"].unique())

# read data on streets
str_all = gpd.read_file(path+'all_streets_ma.shp')
str_all = str_all.to_crs(epsg=26986)

# read data with towns borders
cty_brd = gpd.read_file(path+'counties_ma.shp')
twns = gpd.read_file(path+'towns_ma.shp')
twns = gpd.sjoin(twns, cty_brd, how="inner", op="within")
twns["cty"] = twns['name_right'].str.replace(' County', '')
twns["cty_towns"] = twns["cty"].str.upper() + ', ' + twns["name_left"].str.upper()
# rename Manchester
twns.loc[twns['cty_towns']=='ESSEX, MANCHESTER-BY-THE-SEA','cty_towns'] = 'ESSEX, MANCHESTER'

# create an empty directory for storing the resulting files 
#os.system("rm -rf " + path + "check_massgis_validity && mkdir " + path + "check_massgis_validity")
#cty_towns[cty_towns.index('ESSEX, MANCHESTER')] = 'ESSEX, MANCHESTER-BY-THE-SEA'

for cty_town in cty_towns:
    
    print('processing ' + cty_town + '; ' + str(cty_towns.index(cty_town)+1) + \
    ' out of ' + str(len(cty_towns)))
    # town's name and borders
    twn = twns.iloc[np.where(twns['cty_towns']==cty_town)[0][0],0]
    twn_border = twns.loc[twns['cty_towns']==cty_town,['cty','cty_towns','geometry']]
    twn_border = twn_border.to_crs(epsg=26986)
    twn_border['geometry'] = twn_border['geometry'].buffer(20) # add 20m to towns' borders to include streets that are not within towns, but which have addresses within towns
    
    # read all MassGIS addresses for a given town
    cur_addr = gpd.read_file(path + "cty_towns/mgis_" + cty_town + ".shp")
    # drop records with either missing FULL_NUMBER or STREET_NAME
    cur_addr = cur_addr.loc[(~cur_addr['FULL_NUMBE'].isna()) & \
    (~cur_addr['STREET_NAM'].isna())]
    
    # replace empty string in "UNIT" with NA
    cur_addr.loc[cur_addr["UNIT"]=="",'UNIT']=np.nan
    # drop addresses with "-" (usually, MassGIS has same addresses in separate points w/o "-")
    cur_addr = cur_addr.loc[~cur_addr['FULL_NUMBE'].str.contains('-', regex=False)] 
    # drop addresses with any STATUS other than "active" (for Boston status is all None, hence check with if)
    if any(~cur_addr['STATUS'].isna()):
        cur_addr = cur_addr.loc[cur_addr['STATUS'] == 'ACTIVE']
    
    cur_addr = cur_addr.drop(['STATUS'], axis = 1)
    # replace None in UL_PC_NAME with an empty string, ""
    if any(cur_addr['UL_PC_NAME'].isna()):
        cur_addr.replace({'UL_PC_NAME': None}, "", inplace=True)
    
    # create 'point_id' to identify duplicates
    cur_addr["x"] = cur_addr.apply(lambda cur_addr: cur_addr["geometry"].xy[0][0], axis = 1)
    cur_addr["y"] = cur_addr.apply(lambda cur_addr: cur_addr["geometry"].xy[1][0], axis = 1)
    cur_addr['point_id'] = cur_addr.groupby(['x','y']).ngroup()
    cur_addr = cur_addr.drop(['x','y'], axis=1)
    cur_addr = cur_addr.sort_values(by=["point_id","STREET_NAM","FULL_NUMBE","UNIT"])
    
    # separate "suspicious" addresses -- those addresses that are either
    # 1) too few (only 1/2 addresses with a given street name) + NO such street nearby in OSM
    cur_addr["n_addr"] = cur_addr.groupby(['UL_STREET_','UL_TOWN','UL_PC_NAME'])['MASTER_ADD'].transform('count')
    str_twn = gpd.sjoin(str_all, twn_border, how = 'inner', op='intersects')
    names_str_twn = list(set(str_twn['name']))
    str_twn = str_twn.to_crs(cur_addr.crs) # convert streets' CRS to that of MAD
    
    # 2) too far from the rest of addresses with same street name (more than, say, 2 median distances to nearest neighbor)
    # ... aggregate points with same street names into a cluster
    addr_clust = pd.DataFrame(cur_addr.groupby(['UL_STREET_','UL_TOWN','UL_PC_NAME'])['geometry'].apply(lambda x : x.unary_union)).reset_index()
    addr_clust.rename(index=str, columns={'geometry':'cluster'}, inplace=True)
    cur_addr = cur_addr.merge(addr_clust, on = ['UL_STREET_','UL_TOWN','UL_PC_NAME'], how = 'left')
    # ... subtract a given point from a cluster
    if (cty_town!="SUFFOLK, BOSTON"):
        cur_addr['cluster'] = cur_addr[['cluster','geometry']].apply(lambda x: x['cluster'].difference(x['geometry']), axis=1)
        # ... find the distance from a given point to the remaining points in cluster (minimal distance is returned)
        cur_addr['d_to_oth'] = cur_addr[['cluster','geometry']].apply(lambda x: x['geometry'].distance(x['cluster']), axis=1)
    else: # Data for Boston is too large, need to process it piece-wise
        N = cur_addr.shape[0]
        step = N//200
        i = 0
        start = 0
        while start < N:
            print("Processing ", str(i+1), " out of ~" + str(200))
            start = i*step
            end = min(N, (i+1)*step) - 1
            i += 1
            cur_addr.loc[start:end, 'cluster'] = cur_addr.loc[start:end,['cluster','geometry']].apply(lambda x: x['cluster'].difference(x['geometry']), axis=1)
            # ... find the distance from a given point to the remaining points in cluster (minimal distance is returned)
            cur_addr.loc[start:end,'d_to_oth'] = cur_addr.loc[start:end,['cluster','geometry']].apply(lambda x: x['geometry'].distance(x['cluster']), axis=1)
            cur_addr.loc[start:end, 'cluster'] = None # to avoid memory overload
    
    cur_addr = cur_addr.drop(['cluster'], axis = 1) # drop the 'cluster' variable
    # ... compute median distance from a point to its nearest neighbor in a cluster
    cur_addr['med_d_to_oth'] = cur_addr.groupby(['UL_STREET_','UL_TOWN','UL_PC_NAME'])['d_to_oth'].transform(lambda x : x.median())
    
    mgis_to_check = cur_addr.loc[((cur_addr['d_to_oth'] > 100) &\
        (cur_addr['d_to_oth'] > 10*cur_addr['med_d_to_oth'])) |\
        (   (cur_addr["n_addr"]<4) & ~(cur_addr["UL_STREET_"].isin(names_str_twn))),\
        ['MASTER_ADD', 'FULL_NUMBE', 'POINT_TYPE', 'UL_STREET_',\
        'UL_UNIT', 'UL_BUILDIN', 'UL_COMMUNI', 'UL_TOWN', 'UL_PC_NAME',\
        'geometry', 'point_id', 'd_to_oth','med_d_to_oth','n_addr']]
    
    # Refine the mgis_to_check file by checking the names of nearby streets
    # Remove addresses from mgis_to_check if a street within 100m has the same name as in MAD address
    mad_check_buffer = mgis_to_check.copy()
    mad_check_buffer['geometry'] = mad_check_buffer['geometry'].buffer(100)
    str_twn = str_twn.drop(['index_right'], axis = 1)
    join1 = gpd.sjoin(mad_check_buffer, str_twn, how='left', op='intersects')
    join1['same_name'] = (join1['name']==join1['UL_STREET_'])
    join1 = join1.loc[join1['same_name']]
    mgis_to_check = mgis_to_check.loc[~mgis_to_check["MASTER_ADD"].isin(join1['MASTER_ADD'])]
    
    # Save the resulting file
    if len(mgis_to_check)>0:
        mgis_to_check.to_file(path+'check_massgis_validity/'+cty_town+'_check_MassGIS_validity.shp')

