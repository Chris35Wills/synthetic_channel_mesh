###########
## Run s03_restructure_data.py 

import sys
import s03_restructure_data as rd

path="../test_outputs"
glob_path="*densified_path_*_normals_avg_*m.csv"

rd.restructure(path, glob_path=glob_path)
