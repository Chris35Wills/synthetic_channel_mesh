###########
## Run s05_get_channel_bank_elevation.py 	

import sys
import glob
import s05_get_channel_bank_elevation as gcbe

path="../test_outputs"
obs_path="../test_data"

normal_files=glob.glob("%s/*REARRANGED_CLIPPED.csv" %path)
land_obs="%s/land_obs_xyz.csv" %obs_path
search_dist=10000 # metres

gcbe.get_bank_elevation_DRIVER(land_obs, normal_files, opath=path, search_dist=search_dist)