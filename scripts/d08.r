
###########
## Run s08_get_minimum_surface_from_stack.r

source("./s08_get_minimum_surface_from_stack.r")

raster_path="../test_outputs/"
raster_glob="channel_point2ras_aggregated_channel_*m.tif"
opath="../test_outputs/"

get_min_from_stack(raster_path, raster_glob, opath)


