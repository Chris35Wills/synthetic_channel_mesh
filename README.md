# Synthetic channel mesh development

The scripts and functions available here facilitate the creation of points normal to the centreline of a channel to its edges as constrained by a land classification mask. These points are then assigned elevations. 

<img src="./figs/mesh_normal_method.png" width="300px" />

The long profile elevation trend of the channel is derived from elevations at the head and mouth of the channel, including also any point observations within the channel. Where the centreline itself has only a known elevation at the head and mouth of its length, the long profile elevation trend is simply linear. The cross-profile is constructed as a second order polynomial, using elevations at either edge (using the nearest observed elevations) and the centreline elevation. The meshing procedure is described in [Williams et al., 2017](http://www.the-cryosphere.net/11/363/2017/tc-11-363-2017.html).

Image of example surface before...

<img src="./figs/aoi_without_synth.png" width="500px" />

and after....

<img src="./figs/aoi_with_synth.png" width="500px" />

## System requirements

The code is written in both Python and R. Your Python environment requires the following packages:

	- xxx
	- xxx
	- xxx
	- xxx

These can be installed using the anaconda installer such as by using:

	anaconda ....

R and requires the following libraries to be installed:

	- xxx
	- xxx
	- xxx
	- xxx
	
## Data Requirements

A series of test data files are available in ./test_data which consist of the following:

* centreline xy coordinates (x,y csv file)
* elevation observations (x,y,z csv file)
* land classification grid (channel (1) / not channel (0) raster)

The coordinates within each file must be relative to the same projection - this is not checked.

## How it works

- Centreline points: For a given area of interest (AOI), the user has a pre-digitised channel centreline
- Mask: The centreline falls spatially within a mask of value 1 and not 0 - for an ocean/land mask, the ocean would be 1 and land would be 0

- 01_densify_paths.py: Step XX is to densify the points along the centreline - the denser the points, the denser the output mesh 
- 02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R: Step XX calculates points normal to each centreline points relative the vector between points x units either side of a given centreline node - the larger unit x, the smoother the profile (linearising the channel profile to a greater extent) 
- 03_restructure_data.py: Step XX is a restructure of the data 
- 04_mask_clip.py: Step XX clips the points and their edges using the mask - limits overflow of channel normal points where channel does not have clearly defined sides 
- 05_get_channel_bank_elevation.py: Step XX  <<< observations outside of channel (i.e. bank elevations)
- 06_channel_parabola_EDGE_ELEVATIONS_piecewise.r <<< seed and mouth elevation + observations inside channel
- 07_point_to_raster.r e.g. 14_pnt2ras_E.r 
- 08_get_minimum_surface_from_stack.r 
- 09_combine_and_interp.r

- combine with other points
- interpolate (give some spline code)

## Implementation

	[x]	01_densify_paths.py 		***		
	[x]	02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R 	***
	[x]	03_restructure_data.py 	***
	[x]	04_mask_clip.r 	***
	[x]	05_get_channel_bank_elevation.py 	***
	[x]	06_channel_parabola_EDGE_ELEVATIONS_piecewise.r *** 
	[x]	07_point_to_raster.r + 14_pnt2ras_N.r 
	[x]	08_get_minimum_surface_from_stack.r
	[x] 09_combine_and_interp.r

[ ] Show examples of how to run code at bottom of each script
	-- R equivalent of pythons if __name__ == "__main__" is:
		if (getOption('run.main', default=TRUE)) {
	 		 main()
		}

	see here: https://stackoverflow.com/questions/21383058/is-ifinteractive-an-r-equivalent-to-the-pythonic-if-name-main

[ ] Show how to call code from a single bash script 
	- see ./scripts/DIY.sh
	- use the exampls from the bottom of each script

[ ] Show how to create a quick look spline overview of the improvements...

Explain:
[ ] - true_centre
[ ] - sensitivity to observations - if anything is available, it will tie to it...

## Example run

### Toy data and code

	example_script.sh
	./GitHub/synthetic_channel_mesh/test_data
	./GitHub/synthetic_channel_mesh/test_outputs

	Mask: 		 godthabsfjord_mask__CROP.tif
	Land Obs: 	 land_obs_xyz__clipped.csv
	Channel Obs: channel_obs_xyz__clipped.csv
	Paths:		 path_*.csv (3 are provided)




