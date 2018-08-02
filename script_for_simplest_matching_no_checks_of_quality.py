# The current script matches OSM buildings in MA with address points from MassGIS
# SOME BUILDINGS HAVE NAMES -- keep those variables as well?
import geopandas as gpd
import pandas as pd
import osmnx as ox
import sys
import os
from titlecase import titlecase
#import matplotlib.pyplot as plt

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
    # TO MODIFY: at the very beginning, remove addresses from MassGIS which are already in OSM...
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

