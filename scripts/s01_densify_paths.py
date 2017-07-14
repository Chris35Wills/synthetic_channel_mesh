#!/usr/bin/env python

"""
Program: *_densify_paths.py

Takes in paths output from the pathfinder algorithm and creates additional points 
along the paths at a user set spacing (desired_spacing << units in m)

Function also not to densify but to just reformat data to x and y columns as is the output of densify_driver()

@author: Chris Williams
@date: 04/03/16
"""

import glob 
import os
import pandas as pd
import Linear_referencing as lrf 

def densify_driver(path, file_glob="path_*.csv", desired_spacing=200, opath=''):
	
	filenames=glob.glob("%s/%s" %(path, file_glob))

	if opath=='': 
		opath=path

	for i in range(len(filenames)):
		print("Densifying %s" %(os.path.basename(filenames[i])))

	# xy_cols - the columns of the clipped file that relates to the coord_x and coord_y cols
		lrf.densify(filenames[i], \
			desired_spacing=desired_spacing, \
			xy_cols=[0, 1], \
			plotting=False, \
			out_dir=opath
		)


def reformat_no_densify(path, file_glob="path_*.csv", x="coord_x", y="coord_y"):
	
	filenames=glob.glob("%s/%s" %(path, file_glob))
	
	for i in range(len(filenames)):
	    print("Reformatting but not densifying %s" %(os.path.basename(filenames[i])))

	    # xy_cols - the columns of the clipped file that relates to the coord_x and coord_y cols
	    data=pd.read_csv(filenames[i])
	    ofile=outfile=("%s/NOT_densified_%s.csv" %(path,os.path.splitext(os.path.basename(filenames[i]))[0]))
	    data_out=data[[x, y]]
	    data_out.to_csv(ofile, header=False, index=False)


if __name__ == '__main__':
	
	print("Path densification script - run from import...\n")
	
	#print("Running with test data now - see ../test_outputs/\n")
	#path="../test_data/"
	#file_glob="path_*.csv"
	#opath="../test_outputs/"
	#densify_driver(path, opath=opath, file_glob=file_glob, desired_spacing=100) # << use this to densify
	##reformat_no_densify(path, file_glob=file_glob, x='x', y='y')   # << use this if no densfying needed
