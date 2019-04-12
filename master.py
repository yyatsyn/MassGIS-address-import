import sys

# MODIFY (!) the next string: it must contain a path to a directory with scripts
# (including this script, "master.py") 
sys.path.append('/home/yury/MEGA/OSM/MA addresses')
import my_functions as myf

# MODIFY (!) the file "set_path.py" -- it must contain a directory where all 
# results and temporary files will be stored
from set_path import path

#=======================================================================
# Next two (2) steps must be run to get and prepare the data for imports
#
# Download and process OSM data (takes around 1.5 hours)
import get_osm_data

# Download and process MassGIS address data (takes around 1 hour)
#import get_massgis_data
import get_massgis_data2

#=======================================================================
# Next three (3) steps should be run for some quality checks, but are not
# mandatory for imports
# The order in which they are run doesn't matter

# Find buildings that are in MassGIS, but not in OSM (takes around 2 hours)
import missing_buildings

# Find MassGIS address points for which street names do not match any
# nearby street in OSM
# Fuzzy match ignores words "Street", "Lane", etc., ignores " ", abbreviations and cases
# Runs around 12 hours
import massgis_osm_mismatches_fuzzy

# Find MassGIS address points for which street names do not match any
# nearby street in OSM
# Case insensitive match ignores cases
# Runs around 12 hours
import massgis_osm_mismatches_case_insensitive

#=======================================================================
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
import massgis_check_validity


# Match MassGIS address point to OSM building
# the result is in 7 shp-files per cty_town, stored in path + 'mgis_osm_import_cty_town'
import match_mgis_addr_to_osm_buildings

### To do next
# 0) combine stacked MAD points into one point using ";" (use intervals for units, but not for housenumbers?)
# 1) sjoin MAD and OSM buildings for "BC"-type MAD points  -- 
# 1.1 -- retain both MAD and OSM addresses and figure out if the OSM one needs modification
# 2) for each non-BC MAD point look for nearest N (10, 20, 50?) OSM points with addresses and check if any of them has a matching address. If not, upload the MAD point to OSM

