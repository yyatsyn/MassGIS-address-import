# MassGIS-address-import

This project aims at creating tools for importing OSM addresses in MA
using MassGIS data.  It intends to allow subsequent imports to get
data newly added in the MassGIS database.

The wiki-page for the import: https://wiki.openstreetmap.org/wiki/Import/Catalogue/MassGIS_Addresses 

MEGA folder with files: https://mega.nz/#F!npsVmAiJ!1N69UR3I7E1PMZEQHDwBPQ 

# Prerequisites

 - postgresql (tested with 9.5)
   - hstore (in contrib in the server source, may need manual installation)
 - postgis (tested with 2.5.0beta1)
 - osm2pgsql (tested with 0.96.0)
 - gdal (including python module)
 - python
   - ? (list python modules)
