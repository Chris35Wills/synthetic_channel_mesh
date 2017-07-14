#!/usr/bin/env bash

echo "Run synthetic mesh process for test AOI"
echo "All driver scripts are located in synthetic_channel_mesh/scripts/"

##    ## Run s01_densify_paths.py
##    python d01.py               
##    echo "d01 run complete"
##    echo "	"
##    echo "	"
##    
##    ## Run s02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R
##    Rscript d02.r                
##    echo "d02 run complete"
##    echo "	"
##    echo "	"
##    
##    ## Run s03_restructure_data.py
##    python d03.py               
##    echo "d03 run complete"
##    echo "	"
##    echo "	"
##    
##    ## Run s04_mask_clip.r
##    Rscript d04.r                
##    echo "d04 run complete"
##    echo "	"
##    echo "	"

## Run s05_get_channel_bank_elevation.py
python d05.py               
echo "d05 run complete"
echo "	"
echo "	"

## Run s06_channel_parabola_EDGE_ELEVATIONS_piecewise.r
Rscript d06.r                
echo "d06 run complete"
echo "	"
echo "	"

## Run s07_point_to_raster.r
Rscript d07.r                
echo "d07 run complete"
echo "	"
echo "	"

## Run s08_get_minimum_surface_from_stack.r
Rscript d08.r                
echo "d08 run complete"
echo "	"
echo "	"

## Run s09_combine_and_interp.r
Rscript d09.r
echo "d09 run complete"
echo "	"
echo "	"



