# This file downloads the OSM data from Geofabrik and splits in into several shp-files
# The resulting shp-files are:
# -- counties_ma.shp
# -- towns_ma.shp
# -- addr_buildings_ma.shp
# -- addr_points_ma.shp
# -- all_buildings_ma.shp
# -- noaddr_buildings_ma.shp
# -- all_streets_ma.shp

import geopandas as gpd
import pandas as pd
import os
import urllib.request
from set_path import path

# Download the OSM-file for MA from Geofabrik
ma_geofabrik = "https://download.geofabrik.de/north-america/us/massachusetts-latest.osm.bz2"
bz2_file = path + "massachusetts-latest.osm.bz2"
urllib.request.urlretrieve(ma_geofabrik, bz2_file) 

# Extract the OSM file
unzip_osm_command = "bzip2 -dk -f " + bz2_file
os.system(unzip_osm_command)

# Filter buildings, address points and streets out of the entire OSM file for MA
filter_buildings_command = "osmfilter " + path + "massachusetts-latest.osm --keep-nodes= --keep-relations='building=' --keep-ways='building=' -o=" + path + "buildings_ma.osm"
filter_addr_command = "osmfilter " + path + "massachusetts-latest.osm --keep-relations='building= and addr:housenumber= and addr:street=' --keep-nodes='addr:housenumber= and addr:street=' --keep-ways='building= and addr:housenumber= and addr:street=' -o=" + path + "addr_ma.osm"
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
cty_to_shp = "ogr2ogr -skipfailures " + path + "temp/ " + path + "cty_ma.osm" + " -lco SHPT=POLYGON"  + " -select " + select_var
os.system("rm -rf " + path + "temp && mkdir " + path + "temp") # make sure the "temp" directory is empty
os.system(cty_to_shp)
cty = gpd.read_file(path+'temp/multipolygons.shp')
cty["full_id"] = "r" + cty["osm_id"]
cty.drop("osm_id", axis = 1, inplace=True)
cty.to_file(path+'counties_ma.shp')
# cty = gpd.read_file(path+'counties_ma.shp')

# ... borders of towns
select_vars = "osm_id,name,addr_county"
twn_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "twn_ma.osm" + " -lco SHPT=POLYGON"  + " -select " + select_vars
os.system("rm -rf " + path + "temp && mkdir " + path + "temp") # make sure the "temp" directory is empty
os.system(twn_to_shp)
twn = gpd.read_file(path+'temp/multipolygons.shp')
twn["full_id"] = "r" + twn["osm_id"]
twn.drop("osm_id", axis = 1, inplace=True)
twn.to_file(path+'towns_ma.shp')
# twn = gpd.read_file(path+'towns_ma.shp')

# ... addresses - buildings
select_vars = "osm_way_id,osm_id,addr_housenumber,addr_street" # ":" is replacede with "_"...
addr_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "addr_ma.osm" + " -lco SHPT=POLYGON -nlt PROMOTE_TO_MULTI"  + " -select " + select_vars
os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
os.system(addr_to_shp)
bld_addr = gpd.read_file(path+'temp/multipolygons.shp')
bld_addr = bld_addr[~pd.isnull(bld_addr["addr_stree"]) & ~pd.isnull(bld_addr["addr_house"])]
bld_addr.loc[~pd.isnull(bld_addr["osm_id"]), "full_id"] = "r" + bld_addr.loc[~pd.isnull(bld_addr["osm_id"]), "osm_id"]
bld_addr.loc[~pd.isnull(bld_addr["osm_way_id"]), "full_id"] = "w" + bld_addr.loc[~pd.isnull(bld_addr["osm_way_id"]), "osm_way_id"]
bld_addr.drop(["osm_way_id","osm_id"], axis = 1, inplace=True)
bld_addr.to_file(path+'addr_buildings_ma.shp')
# bld_addr = gpd.read_file(path+'addr_buildings_ma.shp')

# ... addresses - points
select_vars = "osm_id,addr_housenumber,addr_street" 
addr_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "addr_ma.osm" + " -lco SHPT=POINT"  + " -select " + select_vars
os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
os.system(addr_to_shp)
pnt_addr = gpd.read_file(path+'temp/points.shp')
pnt_addr = pnt_addr[~pd.isnull(pnt_addr["addr_stree"]) & ~pd.isnull(pnt_addr["addr_house"])]
pnt_addr["full_id"] = "n" + pnt_addr["osm_id"]
pnt_addr.drop("osm_id", axis = 1, inplace=True)
pnt_addr.to_file(path+'addr_points_ma.shp')
# pnt_addr = gpd.read_file(path+'addr_points_ma.shp')

# ... all buildings
select_vars = "osm_way_id,osm_id,addr_housenumber,addr_street"
buildings_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "buildings_ma.osm"  + " -lco SHPT=POLYGON -nlt PROMOTE_TO_MULTI"  + " -select " + select_vars
os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
os.system(buildings_to_shp)
bld_all = gpd.read_file(path+'temp/multipolygons.shp')
bld_all.loc[~pd.isnull(bld_all["osm_id"]), "full_id"] = "r" + bld_all.loc[~pd.isnull(bld_all["osm_id"]), "osm_id"]
bld_all.loc[~pd.isnull(bld_all["osm_way_id"]), "full_id"] = "w" + bld_all.loc[~pd.isnull(bld_all["osm_way_id"]), "osm_way_id"]
bld_all.drop(["osm_way_id","osm_id"], axis = 1, inplace=True)
bld_all.to_file(path+'all_buildings_ma.shp')
bld_noaddr = bld_all[pd.isnull(bld_all["addr_stree"]) | pd.isnull(bld_all["addr_house"])]
bld_noaddr.to_file(path+'noaddr_buildings_ma.shp')
# bld_noaddr = gpd.read_file(path+'noaddr_buildings_ma.shp')

# ... network of streets (highway=* name=*)
select_vars = "osm_id,name,highway"
streets_to_shp = "ogr2ogr -skipfailures -f 'ESRI Shapefile' " + path + "temp/ " + path + "streets_ma.osm"  + " -lco SHPT=ARC"  + " -select " + select_vars
os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
os.system(streets_to_shp)
str_all = gpd.read_file(path+'temp/lines.shp')
str_all["full_id"] = "w" + str_all["osm_id"]
str_all.drop("osm_id", axis = 1, inplace = True)
str_all.to_file(path+'all_streets_ma.shp')
# str_all = gpd.read_file(path+'all_streets_ma.shp')
