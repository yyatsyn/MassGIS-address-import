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
mgis_input_file = "/home/yury/MEGA/OSM/MA addresses/MassGIS_AddressPoints_and_Locators/MassGIS_Statewide_Address_Points.gdb"
mgis_output_file = "/home/yury/MEGA/OSM/MA addresses/massgis_addresses_selected_variables.shp"
select_vars = "FULL_NUMBER_STANDARDIZED,STREET_NAME,UNIT,BUILDING_NAME,COMMUNITY_NAME,GEOGRAPHIC_TOWN,POSTCODE,PC_NAME,COUNTY"
convert_gdb_to_shp_command = "ogr2ogr -f 'ESRI Shapefile' '" + mgis_output_file + "' '" + mgis_input_file + "' -select " + select_vars
os.system(convert_gdb_to_shp_command)

# To process MassGIS files on a laptop (8Gb RAM), split the mgis_output_file into smaller files by counties&towns
mgis = gpd.read_file(mgis_output_file)
mgis = mgis[~pd.isnull(mgis["COUNTY"]) & ~pd.isnull(mgis["GEOGRAPHIC"])]
mgis["OBJECTID"] = mgis.index
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
filter_str_command = "osmfilter " + path + "massachusetts-latest.osm --keep='highway= and name=' -o=" + path + "streets_ma.osm"
os.system(filter_buildings_command)
os.system(filter_addr_command)
os.system(filter_cty_command)
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


##### The script is work in progress
##### It runs OK up to this point
##### Below are parts of code used for some trial-and-error matching of MassGIS addresss to OSM buildings 



# Split problematic MassGIS addresses by town-county
mgis = gpd.GeoDataFrame()
for cty in ctys:
    cur = gpd.read_file(path + cty + '_problem_mgis.shp')
    mgis = mgis.append(cur)

mgis["cty_town"] = mgis["COUNTY"] + ", " + mgis["GEOGRAPHIC"]
cty_towns = set(mgis["cty_town"].unique())
cty_towns.discard(None)

os.system("rm -rf " + path + "cty_towns && mkdir " + path + "cty_towns")
for cty_town in cty_towns:
    print("Processing " + cty_town + "...")
    cur = mgis.loc[mgis["cty_town"]==cty_town]
    cur.to_file(path + "cty_towns/" + cty_town + ".shp")








#-----------------------------------------------------------------------
# Step 4: check existing addresses, exclude from MassGIS import those that are already in OSM
gc.collect()
mgis = gpd.read_file(mgis_output_file)



counties = ["Barnstable", 
            "Berkshire", 
            "Bristol", 
            "Dukes", 
            "Essex", 
            "Franklin", 
            "Hampden", 
            "Hampshire", 
            "Middlesex", 
            "Nantucket", 
            "Norfolk", 
            "Plymouth", 
            "Suffolk", 
            "Worcester"]

path = "/home/yury/MEGA/OSM/MA addresses/split_by_counties/"
os.chdir(path)

# Save buildings from counties in files
# Made in a separate loop so that it can be run once and then... 
# ... use the downloaded data when debugging/changing the code below
for county in counties:
    print("Processing " + county + "...")
    place_name = county + " County, Massachusetts"
    addr_points = path + "COUNTY_" + county.upper() + ".shp"
    
    # Load from the existing buildings (using OSMNX)
    osm_buildings = ox.buildings_from_place(place_name)
    # Keep only buildings without addresses
    osm_buildings.rename(columns = {'addr:housenumber' : 'housenumber', 'addr:street' : 'street'}, inplace = True)
    osm_buildings = osm_buildings[(osm_buildings.housenumber.isna()) | (osm_buildings.street.isna())]
    osm_buildings['osmid'] = osm_buildings.index
    osmb = osm_buildings[['osmid','geometry']]
    osmb.to_file('buildings_' + county + ".shp")

#=======================================================================
# Match MassGIS points to OSM buildings using several steps:
# 1) exact matching (point within a building and really close to its centroid)
# 2) buffer matching (points within a N-meter proximity of a unique building)
# 3) unmatched points (stop buffer after, say, 5-7 meters)
#    unmatched can be because of 1) addr for parcels, 2) no buildings, 3) buildings already have addr

# Loop across counties to match points to buildings
for county in counties:
    print("Processing " + county + "...")
    place_name = county + " County, Massachusetts"
    addr_points = path + "COUNTY_" + county.upper() + ".shp"
    
    # Load from the existing buildings (using OSMNX)
    osmb = gpd.read_file('buildings_' + county + ".shp")
    
    # Load MassGIS addresses
    mgis = gpd.read_file(addr_points)
    # Make sure both dataframes use the same projections
    osmb = osmb.to_crs(mgis.crs)
    # Add Centroid of buildings to check if the merged address point truely belongs to the building
    #osmb['centroid'] = osmb.centroid
    
    # Intersect MassGIS points with those buildings that do not have proper addr
    join = gpd.sjoin(mgis, osmb, how="inner", op="intersects")
    # find the distance between address point and building's centroid;
    #join["dist"] = join.apply(lambda join: join["geometry"].distance(join["centroid"]), axis = 1)
    # keep matches only if the address point is close to the centroid (distance between the two is < 1m)
    #join = join[join["dist"] < 1]
    
    # subset unmatched address points
    mgis_unm = mgis[~mgis["OBJECTID"].isin(join["OBJECTID"])]
    
    # --- 2 --- Loop over radii of buffers (br), from 0.1 to 3 m
    br = 0.1
    while br <= 3.01:
        print(county + "...buffer size: " + str(br) + " m")
        # replace points with buffer polygons of radius=br
        mgis_unm["geometry"] = mgis_unm.geometry.buffer(br)
        # sjoin buildings with buffers, add field for the distance (br)
        join_buff = gpd.sjoin(mgis_unm, osmb, how="inner", op="intersects")
        #join_buff["dist"] = br
        # append the resulting dataframe
        join = join.append(join_buff)
        # subset unmatched address points
        mgis_unm = mgis[~mgis["OBJECTID"].isin(join["OBJECTID"])]
        # increase radius of buffers
        br += 0.1
    
    # Check if buffer sjoin resulted in same address for different buildings
    # If yes, remove such points from "join" and update mgis_unm(atched)
    join = join.drop_duplicates("OBJECTID", keep = False)
    mgis_unm = mgis[~mgis["OBJECTID"].isin(join["OBJECTID"])]
    
    # Check the number of different street names within buildings
    # If the number is > 1 then exclude all address point within such buildings from "join"
    # Such points should be added as address points
    streets_in_osmb = pd.DataFrame(join.groupby('osmid')['STREET_NAM'].nunique())
    streets_in_osmb.rename(columns={'STREET_NAM':'nstreets'}, inplace=True)
    streets_in_osmb['osmid'] = streets_in_osmb.index
    streets_in_osmb = streets_in_osmb[streets_in_osmb["nstreets"]>1] # buildings with more than 1 street
    # drop address points that are in buildings with multiple addr:street
    join = join[~join["osmid"].isin(streets_in_osmb["osmid"])]
    mgis_unm = mgis[~mgis["OBJECTID"].isin(join["OBJECTID"])]
    
    # Dealing with units: all points with non-empty "unit" will be added as address points,
    # and, hence, excluded from "join"
    join = join[join["UNIT"]==""]
    mgis_unm = mgis[~mgis["OBJECTID"].isin(join["OBJECTID"])]
    
    # Subset data frames to keep only needed variables
    matched = join[["OBJECTID", "FULL_NUMBE", "STREET_NAM", "POSTCODE", "BUILDING_N", "osmid"]]
    
    # Drop duplicates in matched (to make sure addr:housenumber doesn't have 3,3,4,4 etc.')
    matched = matched.drop_duplicates(['FULL_NUMBE', "STREET_NAM",'osmid'], keep= 'last')
    
    # Sort by housenumbers
    matched = matched.sort_values(by=["FULL_NUMBE"])
    
    # Create addr:housenumber for buildings with multiple address points inside
    addr_hn = pd.DataFrame(matched.groupby('osmid')['FULL_NUMBE'].apply(';'.join))
    addr_hn['osmid'] = addr_hn.index
    addr_hn = addr_hn.rename(index=str, columns={"FULL_NUMBE": "addr_hn"})
    matched = matched.merge(addr_hn, left_on = 'osmid', right_on = 'osmid')
    
    # Update the unmatched dataframe
    mgis_unm = mgis[~mgis["OBJECTID"].isin(matched["OBJECTID"])]
    unmatched = mgis_unm[["OBJECTID", "FULL_NUMBE", "UNIT", "STREET_NAM", "POSTCODE", "BUILDING_N", "geometry"]]
    
    # keep only osmid and address fields (no need to track MassGIS points any further)
    matched = matched[['osmid', 'addr_hn', 'STREET_NAM', 'POSTCODE',"BUILDING_N"]]
    matched = matched.drop_duplicates(['osmid'], keep= 'first')
    
    # create a title case street names for both "matched" and "unmatched" dataframes
    # making use of https://pypi.org/project/titlecase/
    def abbreviations(word, **kwargs):
        if word.upper() in ('AT&T', 'JFK', 'YMCA', 'KSC', 'CCC', 'AES', 'FID', 'USMC'):
            return word.upper()
        exclusions = {"O'ROCK": "O'Rock",
                    "MACARTHUR": "MacArthur",
                    "D'AMBROSIO": "D'Ambrosio",
                    "DEWOLFE": "DeWolfe",
                    "LACLEDE": "LaClede",
                    "MACDONALD": "MacDonald",
                    "VAN DE GRAAFF": "Van de Graaff",
                    "D'ANGELO": "D'Angelo",
                    "O'LOUGHIN": "O'Loughin",
                    "DIBENEDETTO": "DiBenedetto",
                    "DECAROLIS": "DeCarolis",
                    "O'DAY": "O'Day",   
                    "MACLAUGHLIN": "MacLaughlin",  
                    "D'AMICO": "D'Amico",
                    "MACDOUGALD": "MacDougald",
                    "O'NEIL": "O'Neil", 
                    "MACMILLAN": "MacMillan",
                    "O'DUNDEE": "O'Dundee",
                    "O'GRADY": "O'Grady",
                    "O'TOOLE": "O'Toole",
                    "DEWITT": "DeWitt",
                    "DEWALT": "DeWalt"}
        if word.upper() in exclusions.keys():
            return exclusions[word.upper()]
    
    #titlecase("DEWALT'S STREET", callback=abbreviations)
    matched["addr_str"] = matched.apply(lambda matched: titlecase(matched["STREET_NAM"], callback=abbreviations), axis = 1)
    unmatched["addr_str"] = unmatched.apply(lambda unmatched: titlecase(unmatched["STREET_NAM"], callback=abbreviations), axis = 1)
    matched["descr"] = matched.apply(lambda matched: titlecase(matched["BUILDING_N"], callback=abbreviations), axis = 1)
    unmatched["descr"] = unmatched.apply(lambda unmatched: titlecase(unmatched["BUILDING_N"], callback=abbreviations), axis = 1)
    #unmatched["addr_unit"] = unmatched.apply(lambda unmatched: titlecase(unmatched["UNIT"], callback=abbreviations), axis = 1)
    unmatched["addr_unit"] = unmatched["UNIT"]
    
    # Save "matched" as zipped csv file
    matched[['osmid', 'addr_hn', 'addr_str', 'descr', 'POSTCODE']].to_csv('matched_' + county, compression='zip')
    
    #.....................
    # TO MODIFY: before merging stacked points among unmatched addresses, 
    # check if there are already similar addresses in the vicinity (say, a 10 m buffer) of each MassGIS point  
    #...
    # For unmatched points -- merge DUPLICATE points if they have same street and housenumber,
    # combine addr:unit as "A,B,1,13" etc.
    ## create 'point_id' to identify duplicates
    unmatched["x"] = unmatched.apply(lambda unmatched: unmatched["geometry"].xy[0][0], axis = 1)
    unmatched["y"] = unmatched.apply(lambda unmatched: unmatched["geometry"].xy[1][0], axis = 1)
    unmatched['point_id'] = unmatched.groupby(['x','y']).ngroup()
    unmatched = unmatched.drop(['x','y'], axis=1)
    unmatched = unmatched.sort_values(by=["point_id","STREET_NAM","FULL_NUMBE","UNIT"])
    
    # split all unmatched points into those w/o units (they will be imported "as is")...
    # ... and those with units (for them a compound addr:unit will be generated for duplicates)
    unmatched_nounit = unmatched[unmatched["addr_unit"]==''] 
    unmatched_unit = unmatched[unmatched["addr_unit"]!=''] 
    
    # create dataframe with combined units for duplicate points
    addr_un = pd.DataFrame(unmatched_unit.groupby(['point_id','FULL_NUMBE','STREET_NAM'])['addr_unit'].apply(';'.join))
    # create dataframe with combined houses for duplicate points
    addr_noun = pd.DataFrame(unmatched_nounit.groupby(['point_id','STREET_NAM'])['FULL_NUMBE'].apply(';'.join))
    
    # merge unmatched_unit and addr_un
    unmatched_unit = unmatched_unit.drop(["addr_unit"], axis=1)
    unmatched_unit = unmatched_unit.merge(addr_un, left_on=['point_id','FULL_NUMBE','STREET_NAM'], right_on=['point_id','FULL_NUMBE','STREET_NAM'])
    unmatched_nounit = unmatched_nounit.drop(['FULL_NUMBE'], axis=1)
    unmatched_nounit = unmatched_nounit.merge(addr_noun, left_on=['point_id','STREET_NAM'], right_on=['point_id','STREET_NAM'])
    
    # replace "&" with "," -- to make addr:unit look uniformly
    unmatched_unit["addr_unit"] = unmatched_unit["addr_unit"].str.replace("&",";")
    unmatched_nounit['FULL_NUMBE'] = unmatched_nounit['FULL_NUMBE'].str.replace("&",";")
    
    # concat the two unmatched dataframes (w/ and w/o units)
    unmatched = unmatched_unit.append(unmatched_nounit)
    
    # Drop duplicates (keep only 1 instance of each point within ['point_id','FULL_NUMBE','STREET_NAM'])
    unmatched = unmatched.drop_duplicates(['point_id','FULL_NUMBE','STREET_NAM'], keep= 'first')
    
    # Drop OBJECTID, UNIT and point_id as they are no longer needed 
    unmatched = unmatched.drop(["OBJECTID", "UNIT",'point_id'], axis=1)
    unmatched = unmatched.rename(index=str, columns={"FULL_NUMBE": "addr_hn"})
    
    # Save the "unmatched" dataframe as zipped csv files
    unmatched.to_csv('unmatched_' + county, compression='zip')
    unmatched.to_file('unmatched_' + county + ".shp")

