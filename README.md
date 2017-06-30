# Synthetic channel mesh development

The scripts and functions available here facilitate the creation of points normal to the centreline of a channel to its edges as constrained by a land classification mask. These points are then assigned elevations. 

<img src="./figs/mesh_normal_method.png" width="300px" />

The long profile elevation trend of the channel is derived from elevations at the head and mouth of the channel, including also any point observations within the channel. Where the centreline itself has only a known elevation at the head and mouth of its length, the long profile elevation trend is simply linear. The cross-profile is constructed as a second order polynomial, using elevations at either edge (using the nearest observed elevations) and the centreline elevation. The meshing procedure is described in [Williams et al., 2017](http://www.the-cryosphere.net/11/363/2017/tc-11-363-2017.html).

Image of example surface before...

- Point plot
- Surface plot

and after....

- Point plot
- Point plot + synth mesh points
- Surface plot

## System requirements

The code is written in R and requires the following libraries to be installed:

## Data Requirements

A series of test data files are available in ./test_data which consist of the following:

* centreline xy coordinates (x,y csv file)
* elevation observations (x,y,z csv file)
* land classification grid (channel (1) / not channel (0) raster)

The coordinates within each file must be relative to the same projection - this is not checked.

## How it works

- For a given area of interest (AOI), the user has a pre-digitised channel centreline
- The centreline falls spatially within a mask of value 1 and not 0 - for an ocean/land mask, the ocean would be 1 and land would be 0

- Step XX is to densfy the points along the centreline - the denser the points, the denser the output mesh (01_densify_paths.py)
- Step XX calculates points normal to each centrleine points relative the vector between points x units either side of a given centreline node - the larger unit x, the smoother the profile (linearising the channel profile to a greater extent) (02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R)
- Step XX is a restructure of the data (04_restructure_data.py)
- Step XX clips the mask... (09_mask_clip.py)
- Step XX 10_get_channel_bank_elevation.py 	***
- Step XX 12c_channel_parabola_EDGE_ELEVATIONS_piecewise.r
- Step XX 14_point_to_raster.r e.g. 14_pnt2ras_E.r 
- Step XX 14-15c_get_minimum_surface_from_stack.r 

- combine with other points
- interpolate (give some spline code)

## Implementation

## Example code






