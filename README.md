# Synthetic channel mesh routine

The code provided enables the creation of points with assigned elevations to populate channels where only the centreline is mapped. This is achieved by creating points normal to the centreline of a channel up to its banks, the extent of which is constrained by a land classification mask. 

<img src="./figs/mesh_normal_method.png" width="300px" />

The long profile elevation trend of the channel is derived from elevations at the head and mouth of the channel, including also any point observations within the channel. Where the centreline itself has only a known elevation at the head and mouth of its length, the long profile elevation trend is simply linear. The cross-profile is constructed as a second order polynomial, using elevations at either edge (using the nearest observed elevations) and the centreline elevation. The meshing procedure is described in [Williams et al., 2017](http://www.the-cryosphere.net/11/363/2017/tc-11-363-2017.html).

Image of example elevation surface before integration of synthetic channel...

<img src="./figs/aoi_without_synth.png" width="500px" />

...after integration of synthetic channel....

<img src="./figs/aoi_with_synth.png" width="500px" />

...with the channel centrelines overlain...

<img src="./figs/aoi_with_synth_and_paths.png" width="500px" />

## System requirements

The code is written in both Python and R. The python portion of the work was developed using Python 3.5.1. You will need the following Python modules:

	- numpy
	- pandas
	- matplotlib
	- pyproj
	- scipy.spatial

If you are using the [anaconda](https://www.continuum.io/downloads) python distribution, you can set up a new python environment along with these modules using:

```python
conda create -n synth_env python=3.5.1 

source activate synth_env # linux
#activate synth_env # windows

conda install -c anaconda numpy=1.13.1
conda install -c anaconda pandas=0.20.2
conda install -c conda-forge matplotlib=2.0.2
conda install -c conda-forge pyproj=1.9.5.1
conda install -c anaconda scipy=0.18.1
```

The R portion of the code was developed using R version 3.2.2 You will need the following libraries to be installed:

	- raster
	- sp
	- ggplot2
	- FNN
	- foreach
	- pgirmess
	- rgdal
	- maptools
	- tools
	- dplyr

The scripts will automatically check for these and try and install them. This may need you to manual install packages depending on your privileges.

## Data Requirements

To take advantage of this code, you should have available:

- channel coordinates in csv format with columns of: x,y
- position and elevation of observations surrounding the channel in csv format with columns of: x,y,z
- any observations within the channel in csv format with columns of: x,y,z OR knowledge of a location at the head and mouth of the channel you propose to create (discussed later)
- a land classification mask as a .tif , the extent of which all points (above) are located 

A series of test data files are available in ./test_data which consist of the following:

* path_00*.csv     - centreline xy coordinates (x,y csv file)
* land_obs_xyz.csv - elevation observations (x,y,z csv file)
* aoi_mask.tif	   - land classification grid (channel (1) / not channel (0) raster)
* dist_mask.tif    - a distance raster (distance of channel pixels from non-channel pixels)
				   - this can be created later on and is provided here for speed

**The coordinates within each file must all be in the same projection - this is not checked. It is assumed that your coordinates are in metres.**

## How it works

The synthetic mesh routine consists of the following scripts (all stored in `./scripts`):

- 01_densify_paths.py 	
- 02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R 	
- 03_restructure_data.py 
- 04_mask_clip.r 
- 05_get_channel_bank_elevation.py 	
- 06_channel_parabola_EDGE_ELEVATIONS_piecewise.r 
- 07_point_to_raster.r + 14_pnt2ras_N.r 
- 08_get_minimum_surface_from_stack.r
- 09_combine_and_interp.r

To summarize, the scripts implement the following:

- 01_densify_paths.py
	- Creates equally spaced points at a user defined interval based on the centreline input (x,y).

- 02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.r
	- Calculates points normal to each centreline points relative the vector between points, *x* units either side of a given centreline node - the larger unit x, the smoother the profile (linearising the channel profile to a greater extent). 
	- This creates the synthetic mesh points either side of a channel's centreline.

- 03_restructure_data.py
	- Restructure the output of 02*.r, specifying channel sides (required for 04*.py onwards).

- 04_mask_clip.py
	- Clips the points and their edges using the land classification mask - this limits overflow of channel normal points where channel does not have clearly defined sides.

- 05_get_channel_bank_elevation.py
	- Assigns elevations to the synthetic mesh points at the edge of the channel based on the nearest observed points (out of the channel). 

- 06_channel_parabola_EDGE_ELEVATIONS_piecewise.r 
	- Assigns elevations to all points within the synthetic mesh.
	- Elevations at the beginning and end of the channel - the *seed* and *mouth* respectively - can be set either using the nearest neighbour from provided observations within the channel or by being specifically declared.
	- For a given cross-section, any other point observations available within the channel are incorporated with the pre-defined edge elevations to calculate a cross channel parabola from which the centreline elevation is extracted.
	- The centreline elevations that are set - the mouth, seed and any other nodes close to points (for which a threshold value of 1000 m is set assuming your coordinates are in metres) - are then used to assign elevations to all other centreline nodes, providing elevations along the entire centreline. Where only the seed and mouth are known, this will be linear, otherwise it will be piecewise.
	- The routine is sensitive to the presence of any observations provided, assuming they are within the search threshold (default is 1000 m assuming your coordinates are in metres).
	- Parabolas are then calculated for each node using the centreline elevation and prior-set edge elevations, providing elevations across each cross section.

- 07_point_to_raster.r 
	- Grids the mesh nodes to a specific resolution, averaging the elevations where multiple cross-sections overlap such as at meanders. 

- 08_get_minimum_surface_from_stack.r 
	- Where multiple synthetic channels have been created within an AOI, overlaps are likely such as at confluences. 
	- This combines all of the rasterised synthetic channels, taking the minimum elevation at overlaps, thus preserving deeper channels within a system.
	- Returned from this script are raster (GeoTIFF) and point xyz (csv) datasets of the combined synthetic dataset.

- 09_combine_and_interp.r
	- This provides an example of how you can integrate the synthetic points with other observations using simple interpolation to provide a DEM with your synthetic channel - a quick way to compare a surface with and without synthetic intervention.
	- More sophisticated interpolation routines should be experimented with - this is just for a visual representation!

Various helper functions are held within the following files, also in `./scripts`:

- Linear_referencing.py
- nearest_neighbour.py
- util.py
- channel_parabola.r
- parabola_funcs_where_obs_available.r
- plotting.r

## Implementation and example run

At the bottom of each script is an example of how to run the code using the provided test data (see `./test_data`). Running on a script by script basis may be preferable as there are numerous settings that can be altered and you may wish to modify/create new functions to suit your specific needs. 

Once you know the settings you wish to apply, a better way to use the code is to use driver scripts which import the functions - examples of these are also provided in `.\scripts` (with a `d` prefix e.g. `d01.py`). Once setup, these can then be called using a single script such as using bash -- an example script is provided to run the example aoi code: `./scripts/example_aoi_run.sh`.

So to run the example, you can use:

```
python d01.py
Rscript d02.r
...etc
```

or

```
sh example_aoi_run.sh
```

Make sure you are within the `./scripts/` folder when you do this.

One method of implementing this command line approach on windows is through [cygwin](https://www.cygwin.com/).

*If you're on windows, you may need to run `dos2unix test_aoi.sh` first to handle newline characters.*

The input and test data outputs for the example run are located respectively in:
	./GitHub/synthetic_channel_mesh/test_data
	./GitHub/synthetic_channel_mesh/test_outputs

**NB/ The code has been set-up such that all scripts and functions work within this directory structure - if you change these locations or perhaps add the functions to your system path, you will need to alter the path declarations in each script accordingly.**

## Comments and further development

The settings currently defined suited the application for which the code was developed but are not definitive as the synthetic meshing conditions will change per application. Furthermore, elements of the code could be made more efficient/more widely applicable - feel free to make developments as you see fit - this can be done through Github [here](https://github.com/Chris35Wills/synthetic_channel_mesh).







