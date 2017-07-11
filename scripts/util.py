from __future__ import division

import os
import os.path

import numpy as np

def trim_constant_rows_cols(img):
	'''
	Removes any rows or columns of a constant value from a image array
	based on a standard deviation of 0
	'''	
	img_goodrows = img[np.std(img, 1) != 0, :]
	img_no_zeros = img_goodrows[:, np.std(img_goodrows, 0) != 0]
	return img_no_zeros


def check_output_dir(filename):
	'''
	Checks if the designated directory name exists, creating it if it doesn't
	'''
	dirname = os.path.dirname(filename)
	if not os.path.isdir(dirname):	
		print("%s DOESN'T exist...\n" % dirname)
		os.makedirs(dirname) 
		print("...but it does now")


def check_if_empty(file_in):
	size=os.path.getsize(file_in)
	if size == 0:
		print("%s is empty" %file_in)
		empty=True
	else:
		empty=False

	return empty


def live_view_z(img):
	'''
	Enables display of z values in active matplotlib show() window
	Example	implementation:
		plt.imshow(np.log(image))
		plt.gca().format_coord = util.live_view_z(np.log(image))
		plt.show()
	See here for more info:
	http://matplotlib.org/examples/api/image_zcoord.html
	'''
	numrows,numcols = img.shape
	def format_coord(x,y):
		col = int(x+0.5)
		row = int(y+0.5)
		if col>=0 and col<numcols and row>=0 and row<numrows:
			z = img[row,col]
			return 'x=%1.4f, y=%1.4f, z=%5.5f'%(x,y,z)
		else:
			return  'x=%1.4f, y=%1.4f' %(x,y)
	return format_coord


def live_view_z_extent(img, extent):
	'''
	Same as live_view_z but works where extent is set in imshow
	Extent must be passed in as a tuple of xmin, xmax, ymin, ymax
	'''
	xmin, xmax, ymin, ymax = extent		
	def format_coord(x,y):
		col = int(x+0.5)
		row = int(y+0.5)
		if col>=ymin and col<ymax and row>=xmin and row<xmax:
			z = img[row,col]
			return 'x=%1.4f, y=%1.4f, z=%5.5f'%(x,y,z)
		else:
			return  'x=%1.4f, y=%1.4f' %(x,y)
	return format_coord


def get_index_1d(array, select_array):
	"""
	Get the index of values from select_array in array

	Variables:
	array        : 	a 1d numpy array
	select array : 	a 1d numpy array
	
	Return:
	1d array of indicies

	Example:
	a=np.array([1,2,3,4,5,6,6,6,6,7,7,8,9,9])
	b=np.array([6,7,9])
	In  : get_index_1d(a, b)
	Out : (array([ 5,  6,  7,  8,  9, 10, 12, 13], dtype=int64),)
	"""
	idx = np.in1d(array.ravel(), select_array.ravel())
	x_indx = np.where(idx.reshape(array.shape))
	return x_indx


def get_index_2d(array, select_array):
	"""
	Get the index of values from select_array in array

	Variables:
	array        : 	a 2d numpy array
	select array : 	a 2d numpy array

	Return:
	array of indicies in x
	array of indicies in y

	Example:
	a=np.array([1,2,3,4,5,6,6,6,6,7,7,8,9,9],[1,2,3,4,5,6,6,6,6,7,7,8,9,9])
	b=np.array([6,7,9])
	In  : get_index_1d(a, b)
	Out : (array([ 5,  6,  7,  8,  9, 10, 12, 13], dtype=int64),)
	"""
	idx = np.in1d(array.ravel(), select_array.ravel())
	y_indx, x_indx = np.where(idx.reshape(array.shape))
	return y_indx, x_indx