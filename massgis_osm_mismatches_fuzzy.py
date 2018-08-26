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
from set_path import path

# Check the correctness of OSM and imported addresses: 
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
str_all = gpd.read_file(path+'all_streets_ma.shp')
str_idx = str_all["geometry"].sindex
# For each existing addr point create a list of 20 nearest streets (actually, streets with rectangulars nearest to the considered point) 
K = 50 # number of streets to pre-select
k = 10 # number of nearest streets to check
# Create a list of nearest streets for each address point
pnt_addr = gpd.read_file(path+'addr_points_ma.shp')
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
problem_pnt_addr3 = problem_pnt_addr.loc[~problem_pnt_addr["full_id"].isin(join["full_id_left"])]
if len(problem_pnt_addr3) > 0:
    problem_pnt_addr3.to_file(path + 'problem_pnt_addr_fuzzy_match.shp')


# ... existing OSM buildings with addresses
# Create a list of nearest streets for each address point
bld_addr = gpd.read_file(path+'addr_buildings_ma.shp')
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
problem_bld_addr3 = problem_bld_addr.loc[~problem_bld_addr["full_id"].isin(join["full_id_left"])]
if len(problem_bld_addr3) > 0:
    problem_bld_addr3.to_file(path + 'problem_bld_addr_fuzzy_match.shp')


# ... imported MassGIS addresses (process by counties)
ctys = pd.read_csv(path+"town_cty_id.csv")
cty_towns = list(ctys["cty_town"].unique())
ctys = list(ctys["cty"].unique())

os.system("mkdir " + path + "cty_towns/problem_mgis_fuzzy")
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
    problem_mgis3 = problem_mgis.loc[~problem_mgis["MASTER_ADD"].isin(join["MASTER_ADD"])]
    
    # Save the resulting Geopandas DataFrames into files
    if len(problem_mgis3) > 0:
        problem_mgis3.to_file(path + "cty_towns/problem_mgis_fuzzy/" + cty_town + '.shp')

# Find missing streets -- as streets that are in problematic address points of MassGIS, but not in OSM
# Combine all problematic MassGIS points
mgis_prblm = gpd.GeoDataFrame()
for cty_town in cty_towns:
    print("Processing " + cty_town)
    try:
        cur = gpd.read_file(path + "cty_towns/problem_mgis_fuzzy/" + cty_town + '.shp')
        mgis_prblm = mgis_prblm.append(cur)
    except:
        print("There are no problem points in " + cty_town + "!")

misstr = mgis_prblm[["STREET_NAM", "GEOGRAPHIC", "COUNTY"]].drop_duplicates()
misstr.columns = ["name", "town", "county"]
#
twn = gpd.read_file(path+'towns_ma.shp')
twn.loc[(twn["name"]=="Nantucket"),"addr_count"] = "Nantucket County"
twn = twn.loc[~pd.isnull(twn["addr_count"])]
twn["town"] = [x.upper() for x in twn["name"]]
twn["county"] = [x.replace(" County", "").strip().upper() for x in twn["addr_count"]]
twn.drop(["name", "addr_count"], axis=1, inplace=True)
# Intersect streets with towns
osmstr = gpd.sjoin(twn, str_all, how="inner", op="intersects")
# Transform street names to upper case
caller = lambda x: x.upper() if type(x) is str else None
osmstr["name"] = [caller(x) for x in osmstr["name"]]
osmstr = osmstr[["name", "town", "county"]].drop_duplicates()
# Merge problem and OSM streets, select those that are only in MassGIS, but not in OSM 
merged = osmstr.merge(misstr, indicator=True, how='outer')
to_add = merged[merged['_merge'] == 'right_only']
to_add.drop(["_merge"], axis=1, inplace=True)
to_add.to_csv(path + "missing_streets_fuzzy_match.csv", index=False)
