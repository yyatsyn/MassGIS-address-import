"""
This file matches MassGIS addresses to OSM buildings
Stacked points with same addresses up to addr:unit // addr:housenumber are
combined into one point (building) with addr:unit=A;B;etc.

 -> cur_addr (all MAD) = cur_addr_bc + cur_addr_pt
 -> cur_addr_bc = unm_cur_addr_bc + join
 -> join = join_to_bld + join_to_pnt
 -> join_to_bld = join_to_bld + join_to_bld_manual (imported addresses that would create duplicates b/c of missing unit/suburb/housename in MAD)
 -> join_to_bld = join_to_bld + mgis_osm_conflict (matched buildings with OSM addr:street different from that in MAD)
 -> join_to_bld = join_to_bld2 = join_to_bld2 (ready for import!) + dup_join_to_bld (imported addresses that would create duplicates b/c of existing OSM data -- need to be checked manually, similar to join_to_bld_manual)
    See also the corresponding ODG flowchart
"""

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

ctys = pd.read_csv(path+"town_cty_id.csv")
cty_towns = list(ctys["cty_town"].unique())
ctys = list(ctys["cty"].unique())

# read data with all osm-buildings and all address points
#osma_all = gpd.read_file(path + 'addr_points_ma.shp')
osmb_all = gpd.read_file(path + 'all_buildings_ma.shp')

# read data with towns borders
cty_brd = gpd.read_file(path+'counties_ma.shp')
twns = gpd.read_file(path+'towns_ma.shp')
twns = gpd.sjoin(twns, cty_brd, how="inner", op="within")
twns["cty"] = twns['name_right'].str.replace(' County', '')
twns["cty_towns"] = twns["cty"].str.upper() + ', ' + twns["name_left"].str.upper()

# create a directory to store "import-ready" and "need_to_check" files
os.system("rm -rf " + path + "mgis_osm_import_cty_town && mkdir " + path + "mgis_osm_import_cty_town")

for cty_town in cty_towns[16:]:
    print('processing ' + cty_town + '; ' + str(cty_towns.index(cty_town)+1) + \
    ' out of ' + str(len(cty_towns)))
    # town's name and borders
    twn = twns.iloc[np.where(twns['cty_towns']==cty_town)[0][0],0]
    twn_border = twns.loc[twns['cty_towns']==cty_town,['cty','cty_towns','geometry']]
    
    # create separate directory for each cty_town (as there will be several files for each cty_town)
    os.system("mkdir '" + path + "mgis_osm_import_cty_town/" + cty_town + "'")
    
    # subset building and address points to those that are within town's borders
    osmb = gpd.sjoin(osmb_all, twns.loc[twns['cty_towns']==cty_town, ['geometry']], how = 'inner')
    osmb = osmb.drop(['index_right'], axis = 1)
    #osma = gpd.sjoin(osma_all, twns.loc[twns['cty_towns']==cty_town, ['geometry']], how = 'inner')
    #osma = osma.drop(['index_right'], axis = 1)
    
    # Convert osma and osmb projection to that of mgis (distance -- in meters)
    #osma = osma.to_crs(epsg=26986)
    osmb = osmb.to_crs(epsg=26986)
    #osma.crs = from_epsg(26986)
    osmb.crs = from_epsg(26986)
    
    # # create a GeoDataFrame with addresses already present in OSM (represent them as points)
    # osmb_addr = osmb[~osmb['addr_house'].isna() & ~osmb['addr_stree'].isna()]
    # osmb_addr['geometry'] = osmb_addr['geometry'].centroid
    # osma = osma.append(osmb_addr) # osm_addresses
    
    # ## split records with ";"/"," in addr_housenumber and addr_unit into separate records
    # #(take into account possibility of interpolation addresses which are added as lines?
    # # there are few cases like that -- maybe it is easier to modify them manually before the import?..)
    # # ... addr_housenumber with ";" OR ","
    # osma_house = osma['addr_house'].str.split('; |,').apply(pd.Series, 1).stack()
    # osma_house.index = osma_house.index.droplevel(-1) # to line up with df's index
    # osma_house.name = 'addr_house'
    # del osma['addr_house']
    # osma = osma.join(osma_house)
    
    # # ... addr_unit with ";" OR ","
    # # replace addr:unit=1-5 as addr_unit=1;2;3;4;5
    # osma['addr_unit'] = osma['addr_unit'].apply(myf.as_list)
    # osma_unit = osma['addr_unit'].str.split('; |,').apply(pd.Series, 1).stack()
    # osma_unit.index = osma_unit.index.droplevel(-1) # to line up with df's index
    # osma_unit.name = 'addr_unit'
    # del osma['addr_unit']
    # osma = osma.join(osma_unit) # <-- data frame with all addresses (from buildings and points, but NOT other ways)
    
    # # convert text variables to upper case as it is easier to compare with MassGIS data
    # #osma['addr_stree'] = osma['addr_stree'].str.upper()
    # #osma['addr_house'] = osma['addr_house'].str.upper()
    # #osma['addr_hou_1'] = osma['addr_hou_1'].str.upper() # addr:housename
    # osma['addr_unit'] = osma['addr_unit'].str.upper()
    
    # read all MassGIS addresses for a given town
    cur_addr = gpd.read_file(path + "cty_towns/mgis_" + cty_town + ".shp")
    # drop records with either missing FULL_NUMBER or STREET_NAME
    cur_addr = cur_addr.loc[(~cur_addr['FULL_NUMBE'].isna()) & \
    (~cur_addr['STREET_NAM'].isna())]
    #mgis = gpd.read_file(path + "cty_towns/mgis_" + cty_town + ".shp")
    #mgis_dup = mgis.loc[mgis.iloc[:,:-1].duplicated(keep=False)]
    
    # replace empty string in "UNIT" with NA
    cur_addr.loc[cur_addr["UNIT"]=="",'UNIT']=np.nan
    # drop addresses with "-" (usually, MassGIS has same addresses in separate points w/o "-")
    cur_addr = cur_addr.loc[~cur_addr['FULL_NUMBE'].str.contains('-', regex=False)] 
    # drop addresses with any STATUS other than "active"
    # ... (those dropped in-ACTIVE points can be imported later when MassGIS
    # geocodes them precisely or added after manual review)
    if (cty_town!='SUFFOLK, BOSTON'):
        cur_addr = cur_addr.loc[cur_addr['STATUS'] == 'ACTIVE']
    else:
        cur_addr['ADDRESS_ID'] = cur_addr['MASTER_ADD']
    
    cur_addr = cur_addr.drop(['STATUS'], axis = 1)
    # replace None in UL_PC_NAME with an empty string, ""
    cur_addr.replace({'UL_PC_NAME': None}, "", inplace=True)
    
    # stack points that have exactly the same geometry
    # For address points -- merge DUPLICATE points if they have same street and housenumber,
    # combine addr:unit as "A,B,1,13" etc.
    ## create 'point_id' to identify duplicates
    cur_addr["x"] = cur_addr.apply(lambda cur_addr: cur_addr["geometry"].xy[0][0], axis = 1)
    cur_addr["y"] = cur_addr.apply(lambda cur_addr: cur_addr["geometry"].xy[1][0], axis = 1)
    cur_addr['point_id'] = cur_addr.groupby(['x','y']).ngroup()
    cur_addr = cur_addr.drop(['x','y'], axis=1)
    cur_addr = cur_addr.sort_values(by=["point_id","STREET_NAM","FULL_NUMBE","UNIT"])
    
    # find duplicates in mgis
    #cur_addr_dup = cur_addr.loc[cur_addr.iloc[:,:-2].duplicated(keep=False),:]
    
    # separate "suspicious" addresses -- those addresses that are either
    # see the file mass_gis_check_validity.py, use the resulting MANUALLY created lists
    # ..of erroneous MAD addresses to exclude from import
    # -->>> [add code here that excludes confirmed erroneous addresses from cur_addr]
    
    # Split MGIS address points into building centroids (add to buildings, if possible)
    # .. and others (to import as address_points)
    if (cty_town!='SUFFOLK, BOSTON'):
        cur_addr_bc = cur_addr.loc[cur_addr['POINT_TYPE']=='BC']
        cur_addr_pt = cur_addr.loc[cur_addr['POINT_TYPE']!='BC']
    else:
        cur_addr_bc = cur_addr
    
    # sjoin massgis (non-PC points) and osm buildings
    join = gpd.sjoin(cur_addr_bc, osmb, how="inner", op="intersects")
    
    # run sjoin with some buffer radius (br), from 0.1 to 3 m around unmatched points
    unm_cur_addr_bc = cur_addr_bc[~cur_addr_bc["ADDRESS_ID"].isin(join["ADDRESS_ID"])] # unmatched MassGIS addresses
    #unm_osmb = osmb[~osmb["full_id"].isin(join["full_id"])] # OSM building that were not matched with MassGIS addresses
    # Loop over radii of buffers 
    br = 0.1
    
    while (br <= 3.01):
        print(cty_town + "...buffer size: " + str(br) + " m")
        # replace points with buffer polygons of radius=br
        unm_cur_addr_bc["geometry"] = unm_cur_addr_bc.geometry.buffer(br)
        # sjoin buildings with buffers, add field for the distance (br)
        #join_buff = gpd.sjoin(unm_cur_addr_bc, unm_osmb, how="inner", op="intersects")
        join_buff = gpd.sjoin(unm_cur_addr_bc, osmb, how="inner", op="intersects")
        #join_buff["dist"] = br
        # append the resulting dataframe
        join = join.append(join_buff)
        # subset unmatched address points
        unm_cur_addr_bc = cur_addr_bc[~cur_addr_bc["ADDRESS_ID"].isin(join["ADDRESS_ID"])]
        #unm_osmb = osmb[~osmb["full_id"].isin(join["full_id"])]
        # increase radius of buffers
        br += 0.1
    
    join['geometry'] = join['geometry'].centroid
    
    # drop unneeded variables (including housename which is usually the name of an amenity, rather than of the building itself)
    keep_vars = ['FULL_NUMBE', 'UNIT', 'UL_STREET_','UL_COMMUNI', 'UL_TOWN', 'UL_PC_NAME',\
    'UL_BUILDIN', 'POSTCODE', 'ADDRESS_ID','geometry','addr_hou_1', 'addr_city',\
    'addr_subur', 'addr_house', 'addr_stree', 'addr_unit', 'point_id', 'full_id']
    join = join[keep_vars]
    
    ### First -- check what addresses can be added directly on buildings, i.e....
    # number of points (by point_id) and number of street names within a building
    pnts_full_id = pd.DataFrame(join.groupby(['full_id'])['point_id','UL_STREET_'].nunique()).reset_index()
    pnts_full_id = pnts_full_id.loc[(pnts_full_id['point_id']==1) & (pnts_full_id['UL_STREET_']==1)]
    # number of buildings with the same ADDRESS_ID
    blds_addr_id = pd.DataFrame(join.groupby(['ADDRESS_ID'])['full_id'].nunique()).reset_index()
    blds_addr_id = blds_addr_id.loc[blds_addr_id['full_id']==1]
    # select only matches with ADDRESS_ID that are not matched to other building and
    #.. with single point_id and street_name
    join_to_bld = join.loc[(join['full_id'].isin(pnts_full_id['full_id'])) & \
    (join['ADDRESS_ID'].isin(blds_addr_id['ADDRESS_ID']))]
    join_to_pnt = join.loc[~(join['full_id'].isin(pnts_full_id['full_id'])) | \
    ~(join['ADDRESS_ID'].isin(blds_addr_id['ADDRESS_ID']))]
    
    # flag addresses that need addr:unit: those addresses, which w/o units
    # .. will be assigned to several buildings (like 685 BROADWAY in Malden)
    join_to_bld['UNIT'] = join_to_bld['UNIT'].fillna("")
    blds_w_un = pd.DataFrame(join_to_bld.groupby(['FULL_NUMBE','UNIT','UL_STREET_'])['full_id'].nunique()).reset_index()
    blds_wo_un = pd.DataFrame(join_to_bld.groupby(['FULL_NUMBE','UL_STREET_'])['full_id'].nunique()).reset_index()
    # blds_w_un that have more than 1 full_id should be imported manually (usually, they should have UNIT=REAR, but MAD doesn't have it)
    # ... or addr:suburb is needed in this case
    join_to_bld_manual = join_to_bld.merge(blds_w_un.loc[blds_w_un['full_id']>1, ['FULL_NUMBE','UNIT','UL_STREET_']], on = ['FULL_NUMBE','UNIT','UL_STREET_'], how='right')
    join_to_bld = join_to_bld.loc[~join_to_bld['ADDRESS_ID'].isin(join_to_bld_manual['ADDRESS_ID'])]
    # addresses where units are necessary to avoid duplicates
    blds_w_un = blds_w_un.merge(blds_wo_un, on = ['FULL_NUMBE','UL_STREET_'])
    blds_w_un = blds_w_un.loc[(blds_w_un['full_id_x']==1)&(blds_w_un['full_id_y']>1),['FULL_NUMBE','UL_STREET_']]
    blds_w_un.drop_duplicates(keep='first', inplace=True)
    blds_w_un['keep_unit']=1
    join_to_bld = join_to_bld.merge(blds_w_un, how='left')
    join_to_bld['keep_unit'] = join_to_bld['keep_unit'].fillna(0)
    
    # Separate addresses that have different streetnames in OSM and MassGIS
    mgis_osm_conflict = join_to_bld.loc[~(join_to_bld['addr_stree'].isna()) & ~(join_to_bld['addr_stree']==join_to_bld['UL_STREET_'])]
    #mgis_osm_conflict2 = join_to_bld.loc[~(join_to_bld['addr_house'].isna()) & ~(join_to_bld['addr_house']==join_to_bld['FULL_NUMBE'])]
    join_to_bld = join_to_bld.loc[~join_to_bld['ADDRESS_ID'].isin(mgis_osm_conflict['ADDRESS_ID'])]
    
    ### Check how added addresses will look in order to check if any duplicates (for buildings' addresses) will be created
    ### If duplicates are created, don't modify the OSM addresses that are identified as duplicates after the possible import
    ### (add such "duplication-creating" points to the "join_to_bld_manual")
    ### some of them might need addr:suburb to be added
    join_to_bld['imp_city'] = np.where(join_to_bld['addr_city'].isna(),\
        join_to_bld['UL_TOWN'], join_to_bld['addr_city'])
    join_to_bld['imp_street'] = np.where(join_to_bld['addr_stree'].isna(),\
        join_to_bld['UL_STREET_'], join_to_bld['addr_stree'])
    # concat existing addr:housenumber with imported MassGIS FULL_NUMBE
    join_to_bld["addr_house"] = join_to_bld["addr_house"].str.replace(",",";")
    join_to_bld['imp_house'] = np.where(join_to_bld['addr_house'].isna(), \
        join_to_bld['FULL_NUMBE'],\
        join_to_bld['addr_house'].str.cat(join_to_bld['FULL_NUMBE'], sep=';'))
    join_to_bld['imp_house'] = join_to_bld['imp_house'].str.split(';').apply(set) # to keep only unique numbers (to avoid "160;160;162" if the existing and imported numbers are the same)
    join_to_bld['imp_house'] = join_to_bld['imp_house'].apply(";".join)
    # concat existing addr:unit with imported MassGIS UNIT
    join_to_bld["addr_unit"] = join_to_bld["addr_unit"].str.replace(",",";")
    join_to_bld['imp_unit'] = join_to_bld['addr_unit']
    join_to_bld['imp_unit'] = join_to_bld['imp_unit'].fillna("")
    join_to_bld['imp_unit'] = np.where((join_to_bld['keep_unit']==1) & (join_to_bld['imp_unit']==""),\
        join_to_bld['UNIT'], join_to_bld['imp_unit'])
    join_to_bld['imp_zip'] = join_to_bld['POSTCODE']
    # concat existing addr:suburb with imported MassGIS COMMUNITY
    join_to_bld['imp_suburb'] = join_to_bld['addr_subur']
    join_to_bld['imp_suburb'] = np.where(join_to_bld['imp_suburb'].isna(),\
        join_to_bld['UL_COMMUNI'], join_to_bld['imp_suburb'])
    # concat existing addr:housename with imported MassGIS UL_BUILDIN -->>> ???
    # .. or maybe not: UL_BUILDIN is usually what should be in "amenity" tag...
    
    # check for duplicated addresses
    join_to_bld2 = join_to_bld[['full_id','imp_city','imp_house','imp_street','imp_unit','imp_zip','ADDRESS_ID']]
    join_to_bld2_hs = join_to_bld2['imp_house'].str.split(';').apply(pd.Series, 1).stack()
    join_to_bld2_hs.index = join_to_bld2_hs.index.droplevel(-1) # to line up with df's index
    join_to_bld2_hs.name = 'imp_house'
    del join_to_bld2['imp_house']
    join_to_bld2 = join_to_bld2.join(join_to_bld2_hs)
    mad_addr_ids = pd.DataFrame(join_to_bld2.groupby(['full_id'])['ADDRESS_ID'].apply(set))
    mad_addr_ids = mad_addr_ids.reset_index()
    mad_addr_ids['ADDRESS_ID'] = mad_addr_ids['ADDRESS_ID'].apply(tuple)
    join_to_bld2.drop(['ADDRESS_ID'], axis = 1, inplace=True)
    join_to_bld2 = join_to_bld2.merge(mad_addr_ids, on='full_id', how='left')
    
    join_to_bld2.drop_duplicates(keep='first', inplace=True)
    dup_join_to_bld = join_to_bld2.loc[join_to_bld2.duplicated(subset= ['imp_city','imp_house','imp_street','imp_unit','imp_zip'],\
        keep=False),:] #remove records that give same addresses to different buildings
    #join_to_bld2.loc[join_to_bld2['imp_house'].str.contains(';', regex=False)] 
    
    # exclude duplicates -- those in dup_join_to_bld -- for manual check/import
    join_to_bld2 = join_to_bld2.loc[~join_to_bld2['full_id'].isin(dup_join_to_bld['full_id'])]
    
    # replace "" in imp_unit with None, stack imp_unit and imp_house into ;-separated lists
    #join_to_bld2['imp_unit'] = join_to_bld2['imp_unit'].where(join_to_bld2['imp_unit']=="", np.nan)
    imp_units = pd.DataFrame(join_to_bld2.groupby(['full_id'])['imp_unit'].apply(';'.join))
    def nat_sort(a):
        """Input ;-separated list, output -- nat-sorted ;-separated list """
        from natsort import natsorted # for sorting strings that have mixed numbers and letters
        c = list(set(a.split(";"))) # to remove duplicates
        c = natsorted(c, alg=ns.IGNORECASE)
        c = ";".join(c)
        return c
    
    imp_units['imp_unit'] = imp_units['imp_unit'].apply(lambda x: nat_sort(x))
    imp_units = imp_units.reset_index()
    imp_house = pd.DataFrame(join_to_bld2.groupby(['full_id'])['imp_house'].apply(';'.join))
    imp_house['imp_house'] = imp_house['imp_house'].apply(lambda x: nat_sort(x))
    imp_house = imp_house.reset_index()
    join_to_bld2.drop(['imp_house','imp_unit'], axis = 1, inplace=True)
    join_to_bld2 = join_to_bld2.merge(imp_units, on='full_id', how='left')
    join_to_bld2 = join_to_bld2.merge(imp_house, on='full_id', how='left')
    
    # delete duplicates -- list of addresses on building ready for import!
    join_to_bld2.drop_duplicates(keep='first', inplace=True)
    
    # Save the resulting files
    try:
        join_to_bld2.to_csv(path + "mgis_osm_import_cty_town/" + cty_town + "/join_to_bld2.csv", index = None, header=True)
    except:
        pass
    
    try:
        dup_join_to_bld.to_csv(path + "mgis_osm_import_cty_town/" + cty_town + "/dup_join_to_bld.csv", index = None, header=True)
    except:
        pass
    
    try:
        mgis_osm_conflict.to_file(path + "mgis_osm_import_cty_town/" + cty_town + "/mgis_osm_conflict.shp")
    except:
        pass
    
    try:
        join_to_bld_manual.to_file(path + "mgis_osm_import_cty_town/" + cty_town + "/join_to_bld_manual.shp")
    except:
        pass
    
    try:
        join_to_pnt.to_file(path + "mgis_osm_import_cty_town/" + cty_town + "/join_to_pnt.shp")
    except:
        pass
    
    if (cty_town!='SUFFOLK, BOSTON'):
        try:
            cur_addr_pt.to_file(path + "mgis_osm_import_cty_town/" + cty_town + "/cur_addr_pt.shp")
        except:
            pass
    
    try:
        unm_cur_addr_bc.to_file(path + "mgis_osm_import_cty_town/" + cty_town + "/unm_cur_addr_bc.shp")
    except:
        pass

