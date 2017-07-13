#!/bin/bash

echo "This bash script will run everything in one fell swoop..."
echo "assuming of course that your Python and R environments are correctly set-up..."
echo " "

## Running python in bash...
echo "Running Python code..."
echo " "
python << END
print("Hello world from bash")
END

## Running r in bash...
echo "Running R code..."
echo " "
Rscript -e 'cat("hello world")'
#source("./06_channel_parabola_EDGE_ELEVATIONS_piecewise.r")

exit