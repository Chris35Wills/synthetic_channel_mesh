import os
import os.path
import numpy as np

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


