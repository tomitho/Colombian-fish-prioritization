# create grass folder
export INDIR=/data/grass_data
mkdir $INDIR 
cd $INDIR

wget "http://www.spatial-ecology.net/ost4sem/exercise/basic_adv_grass/stream_variables.sh"
rstudio stream_variables.sh &
  
#  Download and unzip a DEM from Earthenv.net:

# wget -O  $INDIR/elevation_1KMmd_GMTEDmd.tif  "https://data.earthenv.org/topography/elevation_1KMmd_GMTEDmd.tif"
# gdalinfo $INDIR/elevation_1KMmd_GMTEDmd.tif       # check data

# Create the GRASS GIS data base and enter GRASS:
grass -f --text -c $INDIR/study_area_elevation.tif $INDIR/grass_location

# load in necessary extensions
g.extension  extension=r.stream.basins
g.extension  extension=r.stream.distance
g.extension  extension=r.stream.order
g.extension  extension=r.hydrodem
g.extension  extension=r.stream.snap
g.extension  extension=v.stream.network 

# load in elevation data
r.in.gdal input=$INDIR/study_area_elevation.tif output=elevation --o
r.info 

# Open GUI and visualize the layers:
g.gui wxpython

# Extract drainage direction and stream network
r.watershed  --h  # see help regarding the options and flags
r.watershed  elevation=elevation drainage=drainage   stream=stream  accumulation=accumulation threshold=90  --o

# Get drainage basins (last downstream segment: -l flag)
r.stream.basins  direction=drainage  stream_rast=stream  basins=subbasins --o

# assign stream orders for the river network
r.stream.order --h
r.stream.order  stream_rast=stream elevation=elevation_cond direction=drainage  accumulation=accumulation   stream_vect=streams_v --o
 
# Write files to disk:
v.out.ogr  input=streams_v  output=$INDIR/stream_network_vector_GRASS_90.gpkg   format=GPKG  type=line  --overwrite

# r.out.gdal  input=str_order  output=$INDIR/str_order.tif  type=Int32  nodata=-9999  --o  -c -m   createopt="COMPRESS=LZW,ZLEVEL=9"
r.out.gdal  input=subbasins     output=$INDIR/subbasin_90.tif    type=Int32  nodata=-9999  --o  -c -m   createopt="COMPRESS=LZW,ZLEVEL=9"
r.out.gdal  input=stream  output=$INDIR/stream/stream_network_GRASS_90.tif  type=Int32  nodata=-9999  --o  -c  -m    createopt="COMPRESS=LZW,ZLEVEL=9"

