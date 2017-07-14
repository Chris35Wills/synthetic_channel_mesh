###########
## Run s06_channel_parabola_EDGE_ELEVATIONS_piecewise.r 

source("./s06_channel_parabola_EDGE_ELEVATIONS_piecewise.r")

input_path="../test_data/"
path="../test_outputs/"

maskF=paste0(input_path, "aoi_mask.tif")
land_obsF=paste0(input_path, "land_obs_xyz.csv")
norm_files=Sys.glob(paste0(path, "*REARRANGED_CLIPPED.csv")) ## path numbers expected to be zero indexed otherwise they won't be sorted

#dist=calc_distance_transform(input_path, maskF) # << you can calculate the distance grid here but it is already calculated for speed - takes a long time for a large mask!
dist=paste0(input_path, "dist_mask.tif")

## Example for a specific channel file with an observation file
channel_obsF=paste0(input_path,"obs_xyz.csv")
channel_parabola_DRIVER(norm_files=norm_files,
                       edge_path=path,
                       mask=maskF,
                       land_obs=land_obsF,
                       obs_path=channel_obsF,
                       dist=dist,
                       true_centre=FALSE,
                       dist_threshold=1000,
                       verbose=0)