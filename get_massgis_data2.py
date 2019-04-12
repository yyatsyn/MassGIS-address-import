#Convert MassGIS gdb file into shp, keep only relevant variables, split the data by towns
# ver2: convers variables to upper/lower cases
import geopandas as gpd
import pandas as pd
import os
import urllib.request
from set_path import path
from titlecase import titlecase


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
                "MACKENZIE": "MacKenzie",
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


# Download the MAD-file for MA from 
mad_url= "https://www.dropbox.com/s/fbofk6oy094pn4m/MassGIS_AddressPoints_and_Locators.zip?dl=1"
mad_file = path + "MassGIS_AddressPoints_and_Locators.zip"
urllib.request.urlretrieve(mad_url, mad_file) 

# Extract the MAD file
unzip_mad_command = "unzip -o " + mad_file + " -d " + path
os.system(unzip_mad_command)

# Get a list of counties and towns -- to split the MA file into smaller ones
mgis_input_file = path + "MassGIS_Statewide_Address_Points.gdb"
ctys_output_file = path + "massgis_cty_towns.shp"
select_vars = "GEOGRAPHIC_TOWN,GEOGRAPHIC_TOWN_ID,COUNTY"
convert_gdb_to_shp_command = "ogr2ogr -f 'ESRI Shapefile' '" + ctys_output_file + "' '" + mgis_input_file + "' -lco ENCODING=UTF-8 -select " + select_vars
os.system(convert_gdb_to_shp_command)
ctys = gpd.read_file(ctys_output_file)
ctys = ctys[~pd.isnull(ctys["COUNTY"]) & ~pd.isnull(ctys["GEOGRAPHIC"])]
#mgis["OBJECTID"] = mgis.index # replace with "MASTER_ADD(RESS_ID)"
ctys["cty_town"] = ctys["COUNTY"] + ", " + ctys["GEOGRAPHIC"]
ctys["GEOGRAPH_1"] = [int(x) for x in ctys["GEOGRAPH_1"]]
ctys = ctys[["GEOGRAPHIC", "GEOGRAPH_1", "COUNTY", "cty_town"]].drop_duplicates()
ctys.columns = ["town", "town_id", "cty", "cty_town"]
ctys.to_csv(path+"town_cty_id.csv", index=False)
cty_towns = list(ctys["cty_town"].unique())
ctys = list(ctys["cty"].unique())

if ctys.__contains__(""):
  ctys.pop(ctys.index(""))

# Convert MAD points from gdb into shp file split by counties and towns
select_vars = "MASTER_ADDRESS_ID,FULL_NUMBER_STANDARDIZED,STREET_NAME,UNIT,BUILDING_NAME,COMMUNITY_NAME,GEOGRAPHIC_TOWN,GEOGRAPHIC_TOWN_ID,POSTCODE,PC_NAME,COUNTY"
os.system("rm -rf " + path + "cty_towns && mkdir " + path + "cty_towns")

for cty in ctys:
    print("Converting " + cty)
    cty_output_file = path+ "massgis_addresses_selected_variables_" + cty + ".shp"
    convert_gdb_to_shp_command = "ogr2ogr -where \"COUNTY='" + cty + "'\" " + "'"+ cty_output_file + "' '" + mgis_input_file + "' -lco ENCODING=UTF-8 -select " + select_vars
    os.system(convert_gdb_to_shp_command)
    
    cty_input_file = cty_output_file
    mgis = gpd.read_file(cty_input_file)
    mgis["cty_town"] = mgis["COUNTY"] + ", " + mgis["GEOGRAPHIC"]
    towns = list(mgis["GEOGRAPHIC"].unique())
    for town in towns:
        cty_town = cty + ", " + town
        print("Processing " + cty_town)
        cur = mgis.loc[mgis["cty_town"]==cty_town]
        # Get TownID (zipfiles in basic and advanced address lists are named by Town IDs)
        town_id = cur["GEOGRAPH_1"].unique()[0]
        town_id = str(str(town_id).zfill(3))
        if town_id!="035": # for all cities, except Boston, add "Status" and "Point Type" variables (Boston doesn't have them)
            # create directory to unzip files
            os.system("rm -rf " + path + "temp && mkdir " + path + "temp")
            url_adv_addr = "http://download.massgis.digital.mass.gov/shapefiles/mad/town_exports/adv_addr/AdvancedAddresses_M" + town_id + ".zip"
            url_addr_pts = "http://download.massgis.digital.mass.gov/shapefiles/mad/town_exports/addr_pts/AddressPts_M" + town_id + ".zip"
            zip_adv_addr = path + "temp/" + "AdvancedAddresses_M" + town_id + ".zip"
            zip_addr_pts = path + "temp/" + "AddressPts_M" + town_id + ".zip"
            xls_adv_addr = path + "temp/" + "AdvancedAddresses_M" + town_id + ".xlsx"
            shp_addr_pts = path + "temp/" + "AddressPts_M" + town_id + ".shp"
            
            urllib.request.urlretrieve(url_adv_addr, zip_adv_addr) 
            urllib.request.urlretrieve(url_addr_pts, zip_addr_pts) 
            
            unzip_adv_command = "unzip -o " + zip_adv_addr + " -d " + path + "temp"
            unzip_pts_command = "unzip -o " + zip_addr_pts + " -d " + path + "temp"
            
            os.system(unzip_adv_command)
            os.system(unzip_pts_command)
            
            #... read xlsx-file with an advanced address list, keep ADDRESS_ID and STATUS,
            adv_adr = pd.read_excel(xls_adv_addr)
            adv_adr = adv_adr[["ADDRESS_ID", "STATUS"]]
            #... read shp-file with basic address points, keep ADDRESS_ID and POINT_TYPE
            bas_pnt = gpd.read_file(shp_addr_pts)
            bas_pnt = bas_pnt[["ADDRESS_ID", "POINT_TYPE"]]
            # merge bas_pnt and adv_adr by ADDRESS_ID
            add_var = adv_adr.merge(bas_pnt, left_on='ADDRESS_ID', right_on='ADDRESS_ID', how="outer")
            # merge add_var and cur, keep only addresses that are in MAD
            cur = cur.merge(add_var, left_on="MASTER_ADD", right_on="ADDRESS_ID", how="left")
        else:
            cur["STATUS"] = None
            cur["POINT_TYPE"] = None
        
        # conver all variables into U/L cases
        cur['UL_STREET_NAM'] = cur.apply(lambda cur: titlecase(cur["STREET_NAM"], callback=abbreviations) if pd.notnull(cur['STREET_NAM']) else "", axis = 1)
        cur['UL_UNIT'] = cur.apply(lambda cur: titlecase(cur["UNIT"], callback=abbreviations) if pd.notnull(cur['UNIT']) else "", axis = 1)
        cur['UL_BUILDING_NAME'] = cur.apply(lambda cur: titlecase(cur["BUILDING_N"], callback=abbreviations) if pd.notnull(cur['BUILDING_N']) else "", axis = 1)
        cur['UL_COMMUNITY'] = cur.apply(lambda cur: titlecase(cur["COMMUNITY_"], callback=abbreviations) if pd.notnull(cur['COMMUNITY_']) else "", axis = 1)
        cur['UL_TOWN'] = cur.apply(lambda cur: titlecase(cur["GEOGRAPHIC"], callback=abbreviations) if pd.notnull(cur['GEOGRAPHIC']) else "", axis = 1)
        cur['UL_PC_NAME'] = cur.apply(lambda cur: titlecase(cur["PC_NAME"], callback=abbreviations) if pd.notnull(cur['PC_NAME']) else "", axis = 1)
        #cur.drop(['STREET_NAM', 'UNIT', 'BUILDING_N', 'COMMUNITY_', 'GEOGRAPHIC', 'PC_NAME'], axis=1, inplace=True)
        cur.to_file(path + "cty_towns/mgis_" + cty_town + ".shp")

mgis = None # to free memory
