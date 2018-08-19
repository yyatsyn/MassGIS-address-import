# The current script matches OSM buildings in MA with address points from MassGIS
import geopandas as gpd
import pandas as pd
import sys
import re
import difflib
import os
import gc
from titlecase import titlecase
import urllib.request

path = "/home/yury/Desktop/"

#-----------------------------------------------------------------------
# Step 1:  Convert MassGIS gdb file into shp, keep only relevant variables

# Download the MAD-file for MA from 
mad_url= "https://www.dropbox.com/s/fbofk6oy094pn4m/MassGIS_AddressPoints_and_Locators.zip?dl=1"
mad_file = path + "MassGIS_AddressPoints_and_Locators.zip"
urllib.request.urlretrieve(mad_url, mad_file) 

# Extract the MAD file
unzip_mad_command = "unzip -o " + mad_file + " -d '/home/yury/MEGA/OSM/MA addresses/MassGIS_AddressPoints_and_Locators/'"
os.system(unzip_mad_command)

# Convert MAD points from gdb into shp file
mgis_input_file = "/home/yury/MEGA/OSM/MA addresses/MassGIS_AddressPoints_and_Locators/MassGIS_Statewide_Address_Points.gdb"
mgis_output_file = "/home/yury/MEGA/OSM/MA addresses/massgis_addresses_selected_variables.shp"
select_vars = "MASTER_ADDRESS_ID,FULL_NUMBER_STANDARDIZED,STREET_NAME,UNIT,BUILDING_NAME,COMMUNITY_NAME,GEOGRAPHIC_TOWN,POSTCODE,PC_NAME,COUNTY"
convert_gdb_to_shp_command = "ogr2ogr -f 'ESRI Shapefile' '" + mgis_output_file + "' '" + mgis_input_file + "' -lco ENCODING=UTF-8 -select " + select_vars
os.system(convert_gdb_to_shp_command)

# To process MassGIS files on a laptop (8Gb RAM), split the mgis_output_file into smaller files by counties&towns
mgis = gpd.read_file(mgis_output_file)
mgis = mgis[~pd.isnull(mgis["COUNTY"]) & ~pd.isnull(mgis["GEOGRAPHIC"])]
#mgis["OBJECTID"] = mgis.index # replace with "MASTER_ADD(RESS_ID)"
mgis["cty_town"] = mgis["COUNTY"] + ", " + mgis["GEOGRAPHIC"]
cty_towns = list(mgis["cty_town"].unique())

os.system("rm -rf " + path + "cty_towns && mkdir " + path + "cty_towns")
for cty_town in cty_towns:
    print("Processing " + cty_town + "...")
    cur = mgis.loc[mgis["cty_town"]==cty_town]
    cur.to_file(path + "cty_towns/mgis_" + cty_town + ".shp")

mgis = None

#-----------------------------------------------------------------------
# Step 2:  download the up-to-date OSM data for MA
# 1) split MassGIS addresses by counties (to save on RAM later)
# 2) extract layers with buildings and address-tagged points (to avoid creating duplicates later)
### /usr/share/gdal/2.2/osmconf.ini need modifications to process addr:*** fields and not to include them in other_tags 

# Download the OSM-file for MA from 
ma_geofabrik = "https://download.geofabrik.de/north-america/us/massachusetts-latest.osm.bz2"
bz2_file = path + "massachusetts-latest.osm.bz2"
urllib.request.urlretrieve(ma_geofabrik, bz2_file) 

# Extract the OSM file
unzip_osm_command = "bzip2 -dk -f " + bz2_file
os.system(unzip_osm_command)

# Filter buildings, address points and streets out of the entire OSM file for MA
path = "/home/yury/Desktop/"
filter_buildings_command = "osmfilter " + path + "massachusetts-latest.osm --keep-nodes= --keep-relations= --keep-ways='building=' -o=" + path + "buildings_ma.osm"
filter_addr_command = "osmfilter " + path + "massachusetts-latest.osm --keep-relations= --keep-nodes='addr:housenumber= and addr:street=' --keep-ways='building= and addr:housenumber= and addr:street=' -o=" + path + "addr_ma.osm"
filter_cty_command = "osmfilter " + path + "massachusetts-latest.osm --keep-relations='border_type=county' --keep-nodes-ways= -o=" + path + "cty_ma.osm"
filter_twn_command = "osmfilter " + path + "massachusetts-latest.osm --keep-relations='admin_level=8' --keep-nodes-ways= -o=" + path + "twn_ma.osm"
filter_str_command = "osmfilter " + path + "massachusetts-latest.osm --keep='highway= and name=' -o=" + path + "streets_ma.osm"
os.system(filter_buildings_command)
os.system(filter_addr_command)
os.system(filter_cty_command)
os.system(filter_twn_command)
os.system(filter_str_command)


# Convert filtered osm-files into shp-files,
# filter points with addr_street = None OR addr_housenumber = None (are filtered in as dependencies, i.e. -- entrance of a building with a proper address, etc.)
# ... borders of counties
select_var = "osm_id,name"
cty_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "cty_ma.osm" + " -lco SHPT=POLYGON"  + " -select " + select_var
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp") # make sure the "temp" directory is empty
os.system(cty_to_shp)
cty = gpd.read_file(path+'temp/multipolygons.shp')
cty.to_file(path+'counties_ma.shp')
# cty = gpd.read_file(path+'counties_ma.shp')

# ... borders of towns
select_var = "osm_id,name,addr_county"
twn_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "twn_ma.osm" + " -lco SHPT=POLYGON"  + " -select " + select_var
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp") # make sure the "temp" directory is empty
os.system(twn_to_shp)
cty = gpd.read_file(path+'temp/multipolygons.shp')
cty.to_file(path+'towns_ma.shp')
# cty = gpd.read_file(path+'towns_ma.shp')

# ... addresses - buildings
select_vars = "osm_way_id,addr_housenumber,addr_street" # ":" is replacede with "_"...
addr_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "addr_ma.osm" + " -lco SHPT=POLYGON"  + " -select " + select_vars
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp")
os.system(addr_to_shp)
bld_addr = gpd.read_file(path+'temp/multipolygons.shp')
bld_addr = bld_addr[~pd.isnull(bld_addr["addr_stree"]) & ~pd.isnull(bld_addr["addr_house"])]
bld_addr.to_file(path+'addr_buildings_ma.shp')
# bld_addr = gpd.read_file(path+'addr_buildings_ma.shp')
# ... addresses - points
select_vars = "osm_id,addr_housenumber,addr_street" 
addr_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "addr_ma.osm" + " -lco SHPT=POINT"  + " -select " + select_vars
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp")
os.system(addr_to_shp)
pnt_addr = gpd.read_file(path+'temp/points.shp')
pnt_addr = pnt_addr[~pd.isnull(pnt_addr["addr_stree"]) & ~pd.isnull(pnt_addr["addr_house"])]
pnt_addr.to_file(path+'addr_points_ma.shp')
# pnt_addr = gpd.read_file(path+'addr_points_ma.shp')

# ... all buildings
select_vars = "osm_way_id,addr_housenumber,addr_street"
buildings_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "buildings_ma.osm"  + " -lco SHPT=POLYGON"  + " -select " + select_vars
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp")
os.system(buildings_to_shp)
bld_all = gpd.read_file(path+'temp/multipolygons.shp')
bld_all.to_file(path+'all_buildings_ma.shp')
bld_noaddr = bld_all[pd.isnull(bld_all["addr_stree"]) | pd.isnull(bld_all["addr_house"])]
bld_noaddr.to_file(path+'noaddr_buildings_ma.shp')
# bld_noaddr = gpd.read_file(path+'noaddr_buildings_ma.shp')

# ... network of streets (highway=* name=*)
select_vars = "osm_id,name,highway"
streets_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "streets_ma.osm"  + " -lco SHPT=ARC"  + " -select " + select_vars
os.system("rm -rf /home/yury/Desktop/temp && mkdir /home/yury/Desktop/temp")
os.system(streets_to_shp)
str_all = gpd.read_file(path+'temp/lines.shp')
str_all.to_file(path+'all_streets_ma.shp')
# str_all = gpd.read_file(path+'all_streets_ma.shp')

#-----------------------------------------------------------------------
# Step 3: check the correctness of OSM and imported addresses: 
# -- does the building/address point has a street with the same name as in addr:street in its neighborhood?

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


# ... existing OSM address points
# For each address point we need to find K nearest streets (say, 3-4)
# To avoid computing all possible pair-wise distances between addresses and streets (ways) make use of R-tree spatial index
str_idx = str_all["geometry"].sindex
# For each existing addr point create a list of 20 nearest streets (actually, streets with rectangulars nearest to the considered point) 
K = 50 # number of streets to pre-select
k = 10 # number of nearest streets to check
# Create a list of nearest streets for each address point
pnt_addr = pnt_addr.reset_index(drop=True)
pnt_addr_nearest_str = [list(str_idx.nearest((pnt_addr.loc[i,"geometry"].x, pnt_addr.loc[i,"geometry"].y), K, objects='raw')) for i in range(len(pnt_addr))]
# Compute the actual distance from the address point to pre-selected 20 streets 
problem_addr = []
for i in range(len(pnt_addr)):
    print("Processing address " + str(i) + " out of " + str(len(pnt_addr)))
    # find actual distance from the address point "i" to K nearest streets (aacording to R-tree spatial index)
    dist = []
    for j in range(K):
        dist.append(pnt_addr.loc[i,"geometry"].distance(str_all["geometry"][pnt_addr_nearest_str[i][j]]))
    
    # sort actual distances (from smallest to largest)
    srtd = dist[:]
    srtd.sort()
    # select k nearest streets, convert their names to upper case
    c = [pnt_addr_nearest_str[i][l] for l in [dist.index(srtd[i]) for i in range(k)]]
    nearest_streets_names = set(str_all.loc[c,"name"])
    nearest_streets_names.discard(None)
    
    # Expand abbreviations, remove words like "Street", etc, convert to upper case
    cur_street = remove_simple_mismatches(pnt_addr.loc[i,"addr_stree"], to_remove, abbrevs)
    nearest_streets = [remove_simple_mismatches(x, to_remove, abbrevs) for x in nearest_streets_names]
    
    # If names of the street in addr:street does not match k nearest streets, then store the "i" index
    match_score = [difflib.SequenceMatcher(None, cur_street, x).ratio() for x in nearest_streets]
    if max(match_score) < 0.9:
        problem_addr.append(i)

problem_pnt_addr = pnt_addr.loc[problem_addr]

# Refine the search results:
# rtree sometimes omits nearby streets; use buffer for pre-selected problematic points
problem_pnt_addr2 = problem_pnt_addr.copy()
problem_pnt_addr2["geometry"] = problem_pnt_addr2["geometry"].buffer(.005)
join = gpd.sjoin(problem_pnt_addr2, str_all, how="inner", op="intersects")
join = join.reset_index(drop = True)
join["addr_stree"] = [remove_simple_mismatches(join["addr_stree"].loc[i], to_remove, abbrevs) for i in range(len(join))]
join["name"] = [remove_simple_mismatches(join["name"].loc[i], to_remove, abbrevs) for i in range(len(join))]
join["match_score"] = [difflib.SequenceMatcher(None, join["name"].loc[i], join["addr_stree"].loc[i]).ratio() for i in range(len(join))]
join = join.loc[join["match_score"] >= 0.8]
problem_pnt_addr3 = problem_pnt_addr.loc[~problem_pnt_addr["osm_id"].isin(join["osm_id_left"])]
if len(problem_pnt_addr3) > 0:
    problem_pnt_addr3.to_file(path + 'problem_pnt_addr.shp')


# ... existing OSM buildings with addresses
# Create a list of nearest streets for each address point
bld_addr = bld_addr.reset_index(drop=True)
bld_addr_nearest_str = [list(str_idx.nearest((bld_addr.loc[i,"geometry"].bounds), K, objects='raw')) for i in range(len(bld_addr))]
# Compute the actual distance from the address point to pre-selected 20 streets 
problem_addr = []
for i in range(len(bld_addr)):
    print("Processing address " + str(i) + " out of " + str(len(bld_addr)))
    # find actual distance from the address point "i" to K nearest streets (aacording to R-tree spatial index)
    dist = []
    for j in range(K):
        dist.append(bld_addr.loc[i,"geometry"].distance(str_all["geometry"][bld_addr_nearest_str[i][j]]))
    
    # sort actual distances (from smallest to largest)
    srtd = dist[:]
    srtd.sort()
    # select k nearest streets, convert their names to upper case
    c = [bld_addr_nearest_str[i][l] for l in [dist.index(srtd[i]) for i in range(k)]]
    nearest_streets_names = set(str_all.loc[c,"name"])
    nearest_streets_names.discard(None)
    
    # Expand abbreviations, remove words like "Street", etc, convert to upper case
    cur_street = remove_simple_mismatches(bld_addr.loc[i,"addr_stree"], to_remove, abbrevs)
    nearest_streets = [remove_simple_mismatches(x, to_remove, abbrevs) for x in nearest_streets_names]
    
    # If names of the street in addr:street does not match k nearest streets, then store the "i" index
    match_score = [difflib.SequenceMatcher(None, cur_street, x).ratio() for x in nearest_streets]
    if max(match_score) < 0.9:
        problem_addr.append(i)

problem_bld_addr = bld_addr.loc[problem_addr]

# Refine the search results:
# rtree sometimes omits nearby streets; use buffer for pre-selected problematic points
problem_bld_addr2 = problem_bld_addr.copy()
problem_bld_addr2["geometry"] = problem_bld_addr2["geometry"].buffer(.005)
join = gpd.sjoin(problem_bld_addr2, str_all, how="inner", op="intersects")
join = join.reset_index(drop = True)
join["addr_stree"] = [remove_simple_mismatches(join["addr_stree"].loc[i], to_remove, abbrevs) for i in range(len(join))]
join["name"] = [remove_simple_mismatches(join["name"].loc[i], to_remove, abbrevs) for i in range(len(join))]
join["match_score"] = [difflib.SequenceMatcher(None, join["name"].loc[i], join["addr_stree"].loc[i]).ratio() for i in range(len(join))]
join = join.loc[join["match_score"] >= 0.8]
problem_bld_addr3 = problem_bld_addr.loc[~problem_bld_addr["osm_way_id"].isin(join["osm_way_id"])]
if len(problem_bld_addr3) > 0:
    problem_bld_addr3.to_file(path + 'problem_bld_addr.shp')


# ... imported MassGIS addresses (process by counties)
os.system("mkdir " + path + "cty_towns/problem_mgis")
for cty_town in cty_towns:
    print("Processing " + cty_town + "...")
    mgis = gpd.read_file(path + "cty_towns/mgis_" + cty_town + ".shp")
    mgis = mgis.to_crs(str_all.crs)
    mgis = mgis.reset_index(drop=True)
    # Create a list of nearest streets for each address point
    mgis_nearest_str = [list(str_idx.nearest((mgis.loc[i,"geometry"].x, mgis.loc[i,"geometry"].y), K, objects='raw')) for i in range(len(mgis))]
    # Compute the actual distance from the address point to pre-selected 20 streets 
    problem_addr = []
    for i in range(len(mgis)):
        print(cty_town + ", processing address " + str(i) + " out of " + str(len(mgis)))
        # find actual distance from the address point "i" to K nearest streets (aacording to R-tree spatial index)
        dist = []
        for j in range(K):
            dist.append(mgis.loc[i,"geometry"].distance(str_all["geometry"][mgis_nearest_str[i][j]]))
        
        # sort actual distances (from smallest to largest)
        srtd = dist[:]
        srtd.sort()
        # select k nearest streets, convert their names to upper case
        c = [mgis_nearest_str[i][l] for l in [dist.index(srtd[i]) for i in range(k)]]
        nearest_streets_names = set(str_all.loc[c,"name"])
        nearest_streets_names.discard(None)
        
        if mgis.loc[i,"STREET_NAM"]!=None: # CHECK ADDRESSES with empty street names!!
            # Expand abbreviations, remove words like "Street", etc, convert to upper case
            cur_street = remove_simple_mismatches(mgis.loc[i,"STREET_NAM"], to_remove, abbrevs)
            nearest_streets = [remove_simple_mismatches(x, to_remove, abbrevs) for x in nearest_streets_names]
            
            # If names of the street in addr:street does not match k nearest streets, then store the "i" index
            match_score = [difflib.SequenceMatcher(None, cur_street, x).ratio() for x in nearest_streets]
            if max(match_score) < 0.9:
                problem_addr.append(i)
    
    problem_mgis = mgis.loc[problem_addr]
    
    # Refine the search results:
    # rtree sometimes omits nearby streets; use buffer for pre-selected problematic points
    problem_mgis2 = problem_mgis.copy()
    problem_mgis2["geometry"] = problem_mgis2["geometry"].buffer(.005)
    join = gpd.sjoin(problem_mgis2, str_all, how="inner", op="intersects")
    join = join.reset_index(drop = True)
    join["STREET_NAM"] = [remove_simple_mismatches(join["STREET_NAM"].loc[i], to_remove, abbrevs) for i in range(len(join))]
    join["name"] = [remove_simple_mismatches(join["name"].loc[i], to_remove, abbrevs) for i in range(len(join))]
    join["match_score"] = [difflib.SequenceMatcher(None, join["name"].loc[i], join["STREET_NAM"].loc[i]).ratio() for i in range(len(join))]
    join = join.loc[join["match_score"] >= 0.8]
    problem_mgis3 = problem_mgis.loc[~problem_mgis["OBJECTID"].isin(join["OBJECTID"])]
    
    # Save the resulting Geopandas DataFrames into files
    if len(problem_mgis3) > 0:
        problem_mgis3.to_file(path + "cty_towns/problem_mgis/" + cty_town + '.shp')

