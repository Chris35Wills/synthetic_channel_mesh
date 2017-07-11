#!/usr/bin/env python

"""
Progam: 04_mask_clip.py

Restructures output of 01_densify_paths.py for input to 05_mask_clip.py

@author: Chris Williams
@date: 08/03/16
"""

import sys
import os
import glob
import pandas as pd 
import numpy as np


def restructure(path, glob_path="*densified_path_*_normals_avg_*m.csv"):

	filenames=glob.glob("%s/%s" %(path, glob_path))

	# A pure pandas approach
	for file in filenames:
		
		print("Working on %s" %file)

		data_right=pd.read_csv(file, sep=',')
		data_left=pd.read_csv(file, sep=',')

		#right
		data_right=data_right.drop('x_left',axis=1)
		data_right=data_right.drop('y_left',axis=1)
		data_right['side']=pd.Series(np.array([1]*len(data_right))) # add side id (1 is right)
		data_right=data_right.rename(columns={'x_right': 'x'})
		data_right=data_right.rename(columns={'y_right': 'y'})

		#left
		data_left=data_left.drop('x_right',axis=1)
		data_left=data_left.drop('y_right',axis=1)
		data_left['side']=pd.Series(np.array([2]*len(data_left))) # add side id (2 is left)
		data_left=data_left.rename(columns={'x_left': 'x'})
		data_left=data_left.rename(columns={'y_left': 'y'})

		#concatenate
		data_lr=pd.concat([data_right, data_left])

		#to csv
		out_name=("%s/%s_REARRANGED.csv" %(path, os.path.splitext(os.path.basename(file))[0]))
		data_lr.to_csv(out_name,sep=',', index=False)

	print(" ")
	print("************************")
	print("Re-structuring complete.")
	print("************************")

if __name__ == '__main__':
	
	print("Data restructuring - run from import...\n")
	print("Running with test data now - see ../test_outputs/\n")
	
	path="../test_outputs"
	glob_path="*densified_path_*_normals_avg_*m.csv"
	restructure(path, glob_path=glob_path)

