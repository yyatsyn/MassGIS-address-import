# The current script matches OSM buildings in MA with address points from MassGIS
import geopandas as gpd
import pandas as pd
import os
import urllib.request
from set_path import path

# Read the file with all OSM buildings
bld_osm_all = gpd.read_file(path+'all_buildings_ma.shp')
# Read the file with towns' borders, create "cty_town" variable (e.g. "MIDDLESEX, MALDEN")
twns = gpd.read_file(path+'towns_ma.shp')
twns.loc[(twns["name"]=="Nantucket"),"addr_count"] = "Nantucket County"
twns.loc[(twns["name"]=="Manchester-by-the-Sea"),"name"] = "Manchester"
twns = twns.loc[~pd.isnull(twns["addr_count"])]
twns["town"] = [x.upper() for x in twns["name"]]
twns["county"] = [x.replace(" County", "").strip().upper() for x in twns["addr_count"]]
twns["cty_town"] = twns["county"] + ", " + twns["town"]

# Join spatially towns and buildings (to use later for subset)
# based on location of buildings' centroids
bld_osm_all["way"] = bld_osm_all["geometry"]
bld_osm_all["geometry"] = [x.centroid for x in bld_osm_all["way"]]
bld_osm_twns = gpd.sjoin(bld_osm_all, twns, how="left", op="within")
bld_osm_twns["geometry"] = bld_osm_twns["way"]
bld_osm_twns.drop("way", axis=1, inplace=True)

# Read the file with towns, counties and town_id (from MAD)
twn_id = pd.read_csv(path + "town_cty_id.csv")
# Create a dictionary with town_id and names 
towns_dict = dict(zip(twn_id["town_id"], twn_id["cty_town"]))

## Loop over towns, find MassGIS buildings that are not in OSM
os.system("rm -rf " + path + "absent_buildings && mkdir " + path + "absent_buildings")
os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
count = 1
for town_id, cty_town in towns_dict.items():
    print("Processing building for " + cty_town + ", count " + str(count))
    count += 1
    # Download, unzip and read MassGIS file with buildings footprints for the current town
    url_mgis_town = "http://download.massgis.digital.mass.gov/shapefiles/structures/structures_poly_" + str(town_id) + ".zip"
    bld_mgis_town = path + "temp/structures_poly_" + str(town_id) + ".zip"       
    
    attempts = 10
    while attempts > 0:
        try:
            urllib.request.urlretrieve(url_mgis_town, bld_mgis_town)
        except HTTPError:
            attempts -= 1
            sleep(5)
            continue
        except:
            print(traceback.format_exc())
        break    
    
    unzip_bld_command = "unzip -o " + bld_mgis_town + " -d " + path + "temp"
    os.system(unzip_bld_command)
    bld_mgis_cur = gpd.read_file(path + "temp/structures_poly_" + str(town_id) + ".shp")
    bld_mgis_cur = bld_mgis_cur.to_crs(bld_osm_all.crs)
    
    ####
    # Subset mgis buildings: keep only those whose centroids are in the current cty_town
    bld_mgis_cur["way"] = bld_mgis_cur["geometry"]
    bld_mgis_cur["geometry"] = [x.centroid for x in bld_mgis_cur["way"]]
    bld_mgis_cur1 = gpd.sjoin(bld_mgis_cur, twns.loc[twns["cty_town"]==cty_town,["geometry", "cty_town"]], how="left", op="within")
    bld_mgis_cur2 = bld_mgis_cur1[~pd.isnull(bld_mgis_cur1["cty_town"])]
    bld_mgis_cur2["geometry"] = bld_mgis_cur2["way"]
    bld_mgis_cur2.drop(["way", "index_right"], axis=1, inplace=True)
    ####
    
    # Select OSM buildings for the current town
    bld_osm_cur = bld_osm_twns.loc[bld_osm_twns["cty_town"]==cty_town]
    bld_osm_cur = bld_osm_cur[["full_id_left", "geometry"]]
    # Intersect polygons of MassGIS and OSM building
    join1 = gpd.sjoin(bld_mgis_cur2, bld_osm_cur, how="left", op="intersects")   
    
    # Select only MassGIS buildings that are not in OSM
    join1 = join1[pd.isnull(join1["full_id_left"])]
    miss_bld = bld_mgis_cur2[(bld_mgis_cur2["STRUCT_ID"].isin(join1["STRUCT_ID"]))]
    
    # Refine the results using MAD: some buildings aren't in OSM b/c they have been demolished
    # Keep in miss_bld only those buildings that have address points with 
    # "building centroid" type
    mgiss_addr = gpd.read_file(path + "cty_towns/mgis_" + cty_town + ".shp")
    mgiss_addr = mgiss_addr.to_crs(bld_osm_all.crs)
    try:
        miss_bld2 = gpd.sjoin(miss_bld, mgiss_addr, how="left", op="intersects")
        miss_bld2 = miss_bld2.loc[miss_bld2["POINT_TYPE"]=="BC"]
    except:
        continue
    
    # Save the missing buildings (if any)
    if len(miss_bld2) > 0:
        save_file = path + "absent_buildings/" + cty_town + ".shp"
        miss_bld2.to_file(save_file)


### Merge all new buildings, intersect them with the network of streets 
### to identify the pieces where the streets need some edits
str_all = gpd.read_file(path+'all_streets_ma.shp')
bld_str_int = gpd.GeoDataFrame()
count = 1
for town_id, cty_town in towns_dict.items():
    print("Processing streets intersection for " + cty_town + ", count " + str(count))
    count += 1
    try:
        cur = gpd.read_file(path + "absent_buildings/" + cty_town + ".shp")
        cur_bld_str_int = gpd.sjoin(cur, str_all, how = "left", op = "intersects")
        cur_bld_str_int = cur_bld_str_int.loc[~pd.isnull(cur_bld_str_int["index_right"]),["STRUCT_ID","geometry"]]
        bld_str_int = bld_str_int.append(cur_bld_str_int)
    except:
        continue

bld_str_int.to_file(path + "streets_intersect_new_buildings.shp")

### Merge all new buildings
bld_absent = gpd.GeoDataFrame()
count = 1
for town_id, cty_town in towns_dict.items():
    print("Appending buildings in " + cty_town + ", count " + str(count))
    count += 1
    try:
        cur = gpd.read_file(path + "absent_buildings/" + cty_town + ".shp")
        bld_absent = bld_absent.append(cur)
    except:
        continue

bld_absent.to_file(path + "all_absent_buildings.shp")
