###########
## Run s09_combine_and_interp.r

source("./s09_combine_and_interp.r")

# Example AOI
data_path="../test_data/"
path="../test_outputs/"

pnts1=paste0(data_path, "land_obs_xyz.csv")
pnts2=paste0(path, "min_channel_dem.csv")

#######################
## with synthetic fjord
ras_ofile=paste0(path, "aoi_with_synth.tif")
plot_ofile=paste0(path, "aoi_with_synth.png")

combine_and_interp(pnt_file_paths=c(pnts1, pnts2), 	
				   ofile=ras_ofile,
				   srs='+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

plot_it(ras_path=ras_ofile, ofile=plot_ofile, min=-600, max=1500)

#######################
## without synthetic fjord
ras_ofile=paste0(path, "aoi_without_synth.tif")
plot_ofile=paste0(path, "aoi_without_synth.png")

combine_and_interp(pnt_file_paths=c(pnts1), 	
				   ofile=ras_ofile,
				   srs='+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

plot_it(ras_path=ras_ofile, ofile=plot_ofile, min=-600, max=1500)