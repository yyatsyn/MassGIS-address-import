# MassGIS-address-import

This project aims at creating tools for importing OSM addresses in MA
using MassGIS data.  It intends to allow subsequent imports to get
data newly added in the MassGIS database.

The wiki-page for the import: https://wiki.openstreetmap.org/wiki/Import/Catalogue/MassGIS_Addresses 

MEGA folder with files: https://mega.nz/#F!npsVmAiJ!1N69UR3I7E1PMZEQHDwBPQ 
(some files -- "absent_buildings", "addr_problems_after_fuzzy_matches" and "addr_problems_after_case_insensitive_matches" -- are also uploaded to Dropbox: https://www.dropbox.com/sh/rj0frho8novib13/AAA5CZ9I_lz6M0ybIJad2lNRa?dl=0)

Google Spreadsheet for tracking edits (mismatches in street names, adding absent buildings): https://docs.google.com/spreadsheets/d/1BRMv2iwsg7ZMUiVwtP9JUD5xO8s98ucfVY_1F3DJDfc/edit#gid=1753549676

# Prerequisites

 - postgresql (tested with 9.5)
   - hstore (in contrib in the server source, may need manual installation)
 - postgis (tested with 2.5.0beta1)
 - osm2pgsql (tested with 0.96.0)
 - gdal (including python module)
 - python
   - pandas
   - geopandas
   - titlecase
   - urllib.request
