# MassGIS-address-import

This project aims at creating tools for importing OSM addresses in MA
using MassGIS data.  It intends to allow subsequent imports to get
data newly added in the MassGIS database.

The wiki-page for the import: https://wiki.openstreetmap.org/wiki/Import/Catalogue/MassGIS_Addresses 

MEGA folder with files: https://mega.nz/#F!npsVmAiJ!1N69UR3I7E1PMZEQHDwBPQ 

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
