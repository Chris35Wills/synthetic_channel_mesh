library(akima)
library(raster)
library(rasterVis)

#' pnt_file_paths - each file expected to be a csv (header x,y,z)
combine_points<-function(pnt_file_paths){
	
	xyz=data.frame()
	expected_colnames=c('x','y','z')
	for (f in pnt_file_paths){
		
		f_xyz=read.csv(f)
		
		if (all((colnames(f_xyz)==expected_colnames)==TRUE)==TRUE){ # check cols - x,y,z
			xyz=rbind(xyz,f_xyz)
		} else{
			stop("All input point files must have three columns of the format: 'x', 'y', 'z'")
		}
	}

	return(xyz)

}

#' xyz - a dataframe with colnames of x, y and z
#' srs - spatial reference system provided proj.4 syntax 
#'     - if declared, used to set the reference system of the output raster
#' nx  - x dimension of grid to inteprolate to
#' ny  - y dimension of grid to inteprolate to
interp_surface<-function(xyz, ofile, srs='', nx=400, ny=400){
	
	# Quik interp using the akima package
	s=interp(xyz$x,xyz$y,xyz$z, nx=nx, ny=ny)
	ras=raster(s)
	
	# Project raster
	if (crs!=''){
		crs(ras)=srs
	}

	# Save raster
	writeRaster(ras, ofile)
}

#' Single call to use combine_points() and interp_surface()
combine_and_interp<-function(pnt_file_paths, ofile, srs='', nx=400, ny=400){
	xyz=combine_points(pnt_file_paths)
	interp_surface(xyz, ofile, srs=srs, nx=nx, ny=ny)
}


plot_it<-function(ras_path, ofile='', min='', max=''){
	
	ras=raster(ras_path)

	# If ofile not set, create file in working directory with timestamp
	if (ofile==''){
		ofile=paste0("./ras_plot_", gsub(":", ".", format(Sys.time(), "%X")), ".png")
		}

	if (min==''){min=cellStats(ras, stat='min')}
	if (max==''){max=cellStats(ras, stat='max')}
	#ofile=paste0(path, "aoi_without_synth.png")
	png(ofile)
	p=levelplot(ras, 
				at=seq(min, max, length=15), 
				par.settings = RdBuTheme,
				margin=FALSE)
	print(p)
	dev.off()

}


if (getOption('run.main', default=TRUE)) {
	print("Run from import ... now running code with example code (will fail if earlier scripts in the processing chain have not already been run)")
	
	# Example AOI
	data_path="../test_data/"
	path="../test_outputs/"

	pnts1=paste0(data_path, "land_obs_xyz.csv")
	pnts2=paste0(path, "min_channel_dem.csv")


	## with synthetic fjord
	ras_ofile=paste0(path, "aoi_with_synth.tif")
	plot_ofile=paste0(path, "aoi_with_synth.png")
	
	combine_and_interp(pnt_file_paths=c(pnts1, pnts2), 	
					   ofile=ras_ofile,
					   srs='+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
	
	plot_it(ras_path=ras_ofile, ofile=plot_ofile, min=-600, max=1500)


	## without synthetic fjord
	ras_ofile=paste0(path, "aoi_without_synth.tif")
	plot_ofile=paste0(path, "aoi_without_synth.png")

	combine_and_interp(pnt_file_paths=c(pnts1), 	
					   ofile=ras_ofile,
					   srs='+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

	plot_it(ras_path=ras_ofile, ofile=plot_ofile, min=-600, max=1500)
}
