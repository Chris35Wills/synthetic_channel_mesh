## Run s01_densify_paths.py

import sys
import s01_densify_paths as dp

path="../test_data/"
file_glob="path_*.csv"
opath="../test_outputs/"

dp.densify_driver(path, opath=opath, file_glob=file_glob, desired_spacing=100)




