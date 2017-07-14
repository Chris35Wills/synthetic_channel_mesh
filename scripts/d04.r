###########
## Run s04_mask_clip.r 

source("./s04_mask_clip.r")

path="../test_outputs/"
maskF="../test_data/aoi_mask.tif"
mask=raster(maskF)
glob_path="*REARRANGED.csv"

mask_clipper(path, mask, limit_different_edge_lengths=TRUE, diff_factor=2, glob_path=glob_path)