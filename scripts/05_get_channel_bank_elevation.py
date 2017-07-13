#!/bin/usr/python

"""
Program : get_channel_bank_elevation.py

Assigns an elevation value to each of the edge points of each channel normal line 
created as part of the channel profile mesh program (02....), which is then modified using 
04_restructure_data.py and 05_mask_clip.py. This script should be run after the mask clip 
and before the parabola is created.

Elevations are assigned to the edge points by considering the elevation of the nearest 
neighbour from the bamber DEM (as points) to the channel edge points. This will ensiure that 
when creating the parabola for each normal in the channel mesh (06_channel_parabola.py), the 
edge points will tie in with the exitsing DEM.

1. Reads in the Bamber DEM points (AOI specific) |	 THESE SHOULD BY THIS POINT BE A SINGLE 
2. Reads in the Morlighem points (AOI specific)  |   DATASET (see 07_get_morlighem_points_within aoi.py)

IMPORTANT: The bamber points are clipped by the mask - these "points" are extracted from the bamber bed 
surface interpolation - despite being clipped, there may still be some inetrpolation bleed at the edges 
(the points of which have been influenced by the overall interpolation steps applied (Bamber et al., 2013)
- consequently, these edges may be lower than the GIMP altimetry points that are available at these bank 
edge locations - it may be prudent to use these in future instead of the Bamber DEM (or at least us an 
updated  version of the Bamber DEM and have some idea of the offset between teh GIMP DEM on land and the 
Bamber DEM product)

3. Calculates the nearest neighbour
4. Appends elevations to end points (nothing appended if not an endpoint)

NB/ end points will be the max (and) min distance from the centreline for a given ID

@author : Chris Williams
@date   : 30/03/16
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nearest_neighbour # ~github\roll_your_own_py\src\geospatial
#import georaster         # ~github\python_functions
import util

def test_length(dists,indxs, infinity=float("inf")):
	"""
	Tests to see that all points find a neighbour within the user set neighbourhood distance window
	Fails if neighbours are not found (for which distance values of Inf are returned)
	"""
	A=indxs[dists!=infinity]
	try:
		assert len(A)==len(indxs)
	except AssertionError:
		sys.exit("Neighbours not being found for all edge points - try a larger neighbourhood search distance")


# Read in the Bamber DEM points (AOI specific) |   THESE SHOULD BY THIS POINT BE A SINGLE 
# Read in the Morlighem points (AOI specific)  |   DATASET (see 12_get_morlighem_points_within aoi.py)

def get_bank_elevation_DRIVER(land_obs, normal_files, search_dist=30000, opath=''):
	"""
	search_dist - set this to the max dist you are interested in - points further than this return a dist value of Inf (and can be culled)
	"""

	if opath=='':
		# if output path not set, save where normal files are saved
		opath=os.path.dirname(normal_files[0]) 

	land_obs="C:/Github/synthetic_channel_mesh/test_data/land_obs_xyz.csv"
	elev_xyz=pd.read_csv(land_obs, sep=",")
	elev_xy=elev_xyz[['x','y']]

	elev_xyz=np.asarray(elev_xyz)
	elev_xy=np.asarray(elev_xy)

	# loop through normal files and associate them with nn elevation points 

	for f in normal_files:
		print("File: %s" %f)
		
		norm=pd.read_csv(f) # expected column headers: x y cl_x cl_y  cl_id  cl_dist_r cl_dist_l  path  side  mask_value

		norm_x=norm['x'].values
		norm_y=norm['y'].values
		norm_cl_x=norm['cl_x'].values
		norm_cl_y=norm['cl_y'].values
		norm_cl_id=norm['cl_id'].values
		norm_cl_dist_r=norm['cl_dist_r'].values
		norm_cl_dist_l=norm['cl_dist_l'].values
		norm_path=norm['path'].values
		norm_side=norm['side'].values

		#create empty array to hold edge elevs
		edge_elev=np.asarray([None]*len(norm_side))

		##Get edge points out of dataframe
		#max norm_cl_dist_l for each cl_id and its row index
		norm_left=norm[norm.side==2]
		edge_LEFT=norm.loc[norm_left.groupby('cl_id').cl_dist_l.agg('idxmax')]
		edge_LEFT_indx=edge_LEFT.index

		#max norm_cl_dist_r for each cl_id and its row index
		norm_right=norm[norm.side==1]
		edge_RIGHT=norm.loc[norm_right.groupby('cl_id').cl_dist_r.agg('idxmax')]
		edge_RIGHT_indx=edge_RIGHT.index

		# as numpy arrays ... as I know how to index them :)
		edge_LEFT_array=np.asarray(edge_LEFT)
		edge_RIGHT_array=np.asarray(edge_RIGHT)

		edge_LEFT_xy=edge_LEFT_array[:,0:2]
		edge_RIGHT_xy=edge_RIGHT_array[:,0:2]

		## Get the nearest neighbours from the elev_xy points (bed elev) to the edge points
		print("Calculating nn for the LEFT edge points...")
		dists_RIGHT,indxs_RIGHT = nearest_neighbour.get_nn(elev_xy, edge_RIGHT_xy, nn=1, radius=search_dist)
		print("Calculating nn for the RIGHT edge points...")
		dists_LEFT,indxs_LEFT = nearest_neighbour.get_nn(elev_xy, edge_LEFT_xy, nn=1, radius=search_dist)

		##Keep only points within threshold (all edge points should have a neighbour!!)
		test_length(dists_RIGHT, indxs_RIGHT)
		test_length(dists_LEFT, indxs_LEFT)

		## get the z value associated with the nearest neighbours
		edge_elev_LEFT=elev_xyz[indxs_LEFT,2]
		edge_elev_RIGHT=elev_xyz[indxs_RIGHT,2]
		##reshape the arrays to ensure they can be stacked with the rest of the edge data
		edge_elev_LEFT=edge_elev_LEFT.reshape(len(edge_elev_LEFT),1)
		edge_elev_RIGHT=edge_elev_RIGHT.reshape(len(edge_elev_RIGHT),1)

		## append this elevation as a new column to the row entry for a given edge point
		edge_LEFT_with_ELEV=np.hstack((edge_LEFT_array, edge_elev_LEFT))
		edge_RIGHT_with_ELEV=np.hstack((edge_RIGHT_array, edge_elev_RIGHT))


	#	********* CHECK THIS CATCH *************
		##check arrays cover all cl positions (or at least the same positions)
	#	l_cl_check=np.unique(edge_LEFT_with_ELEV[:,4])
	#	r_cl_check=np.unique(edge_RIGHT_with_ELEV[:,4])
	#	try:
	#		assert np.array_equal(l_cl_check, r_cl_check)
	#	except AssertionError:
	#		sys.exit("***************\n Arrays not equal - edge points either side not dealing with all of the same centreline nodes \n***************")

		##check left and right arrays only deal with left and right points
		l_side_check=np.unique(edge_LEFT_with_ELEV[:,8]) 
		r_side_check=np.unique(edge_RIGHT_with_ELEV[:,8])
	
	#	if np.array_equal(l_side_check, r_side_check):
	#		sys.exit("arrays shouldn't be equal!!!")

		edges_with_ELEV=np.vstack((edge_LEFT_with_ELEV, edge_RIGHT_with_ELEV))

		##write out as a new csv
		#fout=("%s_edge_elevs.csv" %os.path.splitext(f)[0])
		fout=("%s/%s_edge_elevs.csv" %(opath,os.path.basename(os.path.splitext(f)[0])))
		util.check_output_dir(fout)
		print("Writing out edge points and elevations to file...")	

		#np.savetxt(fout, edges_with_ELEV, delimiter=",", header="x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side,bed elev(m)",fmt='%f,%f,%f,%f,%i,%i,%i,%i,%i,%f')

		#######################################
		# The edge files now contain the mask value (in the 10th column)
		# The next step eliminates that from the output

		edges_with_ELEV__no_elev_or_mask_value=edges_with_ELEV[:,0:9] # Ignore mask value column -- final column of this ( edges_with_ELEV__no_elev_or_mask_value[:,8]) is 'side'
																	  # Assumes columns are of format: x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side,mask_value,bed elev(m)
		edges_with_ELEV__elev=edges_with_ELEV[:,10] # Keep only the elevations
													# This assumes columns are of format: x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side,mask_value,bed elev(m)
		
		# convert to dataframe
		edges_with_ELEV__no_elev_or_mask_value=pd.DataFrame(edges_with_ELEV__no_elev_or_mask_value)
		edges_with_ELEV__elev=pd.DataFrame(edges_with_ELEV__elev)

		edges_with_ELEV__no_elev_or_mask_value.columns=['x','y','cl_x','cl_y','cl_id','cl_dist_r','cl_dist_l','path','side']
		edges_with_ELEV__elev.columns=['bed elev(m)']

		edges_with_ELEV__no_mask_value=edges_with_ELEV__no_elev_or_mask_value
		edges_with_ELEV__no_mask_value['bed elev(m)']=edges_with_ELEV__elev	

		edges_with_ELEV__no_mask_value.to_csv(fout, index=False)			
		
		print("Output file: %s\n" %fout)	

#		bed_dem_neighbours_RIGHT=elev_xyz[indxs_RIGHT,:]	
#		#fout_R=("%s_bed_DEM_neighbours_RIGHT.csv" %os.path.splitext(f)[0])
#		fout_R=("%s/%s_bed_DEM_neighbours_RIGHT.csv" %(opath,os.path.basename(os.path.splitext(f)[0])))
#		util.check_output_dir(fout_R)
#		print("Writing out RIGHT neighbours to file...")	
#		np.savetxt(fout_R, bed_dem_neighbours_RIGHT, delimiter=",", header="x,y,z",fmt='%f,%f,%f')
#
#		bed_dem_neighbours_LEFT=elev_xyz[indxs_LEFT,:]
#		#fout_L=("%s_bed_DEM_neighbours_LEFT.csv" %os.path.splitext(f)[0])
#		fout_L=("%s/%s_bed_DEM_neighbours_LEFT.csv" %(opath,os.path.basename(os.path.splitext(f)[0])))
#		util.check_output_dir(fout_L)
#		print("Writing out LEFT neighbours to file...\n")	
#		np.savetxt(fout_L, bed_dem_neighbours_LEFT, delimiter=",", header="x,y,z",fmt='%f,%f,%f')

	print("************************************")	
	print("Completed edge elevation assignment.")
	print("************************************")	

if __name__ == "__main__":

	print("Run from import...\n")
	print("Running with test data now - see ../test_outputs/\n")
	
	search_dist=2000 # metres

	path="../test_outputs"
	obs_path="../test_data"

	normal_files=glob.glob("%s/*REARRANGED_CLIPPED.csv" %path)
	land_obs="%s/land_obs_xyz.csv" %obs_path

	get_bank_elevation_DRIVER(land_obs, normal_files, opath=path)

