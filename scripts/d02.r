###########
## Run s02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.R 	

source("./s02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.r")

in_path="../test_outputs/"
out_path=in_path
file_glob="densified_*.csv"

# a smaller value dist_between_norm_vertices (m) will provided a coarser mesh but will speed things up
vector_from_norms(in_path, out_path, file_glob=file_glob, dist_between_norm_vertices=200, verbose=FALSE)

