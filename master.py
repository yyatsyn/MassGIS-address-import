import sys

# MODIFY (!) the next string: it must contain a path to a directory with scripts
# (including this script, "master.py") 
sys.path.append('/home/yury/MEGA/OSM/MA addresses')

# MODIFY (!) the file "set_path.py" -- it must contain a directory where all 
# results and temporary files will be stored
from set_path import path 

#=======================================================================
# Next two (2) steps must be run to get and prepare the data for imports
#
# Download and process OSM data (takes around 1.5 hours)
import get_osm_data

# Download and process MassGIS address data (takes around 1 hour)
import get_massgis_data

#=======================================================================
# Next three (3) steps should be run for some quality checks, but are not
# mandatory for imports
#
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

