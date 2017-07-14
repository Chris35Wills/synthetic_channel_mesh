###########
## Run s07_point_to_raster.r + 14_pnt2ras_N.r 

source("./s07_point_to_raster.r")

mask_path="../test_data/"
path="../test_outputs/"
glob_extension="*___PIECEWISE.csv"
opath=path
aoi_mask=paste0(mask_path, "aoi_mask.tif")

points_to_raster(path, mask_path, aoi_mask, opath, glob_extension=glob_extension, aggFactors=c(1))
