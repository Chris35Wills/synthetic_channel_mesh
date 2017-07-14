# Program: *_get_minimum_surface_from_stack.r
# Create a new raster just of the minimum values of all single channel rasters, compiled into a stack
# Min surface should elimnate ridges where channels overlap
#
# @author Chris Williams
# @date 27/04/16

if (!require("raster")) install.packages("raster")

check_rasters_for_all_zeros<-function(ras_path, ras_glob_name){
	
	rasList=Sys.glob(paste0(ras_path, ras_glob_name))
	empty_rasters=data.frame(paths_to_ignore=character())

	for (ras in rasList){
		
		a_raster=raster(ras)
		if (minValue(a_raster) == maxValue(a_raster)){  
			empty_rasters=rbind(empty_rasters, ras) 
		}
	}

	if (nrow(empty_rasters)!=0){
	
		cat("Some rasters contain only the same values - these will be ignored in the aggregation...\n")		
	
		empty_raster_df=as.data.frame(empty_rasters)
		rasList_df=as.data.frame(rasList)
		
		colnames(empty_raster_df)<-c('path_to_ignore')
		colnames(rasList_df)<-c('valid_path')

		valid_rasters=rasList[!(rasList_df$valid_path %in% empty_raster_df$path_to_ignore)]
		
		return(valid_rasters)	

	} else {

		return(rasList)

	}
}

get_min_from_stack<-function(raster_path, raster_glob, opath){
	
	raster_path_list=Sys.glob(paste0(raster_path, raster_glob))
	ras_stack=stack(raster_path_list) 

	#define raster projections
	coords<-"+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
	projection(ras_stack) <- coords

	#get raster of minimum values of all cells in the raster stack
	stack_min <- min(ras_stack, na.rm=TRUE)
		
	#save raster
	#ras_out=capture.output(cat(opath, "min_channels_150_m_1200m_every_20_singChannelApprch.tif", sep=""))
	ras_out=capture.output(cat(opath, "min_channel_dem.tif", sep=""))
	writeRaster(stack_min, ras_out, overwrite=TRUE)

	cat("\n\n\nExtracted minimum surface from stack\n")
	cat("No weighting is applied - all surfaces were treated as equal\n")
	cat("I've saved it here:", ras_out, "\n")

	pnts=rasterToPoints(stack_min)
	colnames(pnts)<-c('x','y','z')
	pnts_out=capture.output(cat(opath, "min_channel_dem.csv", sep=""))
	write.csv(pnts, pnts_out, row.names=FALSE)
	cat("I've saved the points here:", pnts_out, "\n")
}

if (getOption('run.main', default=TRUE)) {
	print("Run from import ...")

	#print("Now running code with example code (will fail if earlier scripts in the processing chain have not already been run)")
	#
	#raster_path="../test_outputs/"
	#raster_glob="channel_point2ras_aggregated_channel_*m.tif"
	#opath="../test_outputs/"
	#get_min_from_stack(raster_path, raster_glob, opath)
	#
	#############
	### Plot data
	#library(plot3D)
	#ras=raster("../test_outputs/min_channel_dem.tif")
	#ras=aggregate(ras, fact=5)
	#zData <-round(as.matrix(ras),1)
	#
	#ofile_name=capture.output(cat(opath, "min_channel_3d.png", sep=''))
	#png(ofile_name)
	#persp3D(z=zData, xlab = "easting", bty = "g",
	#  	  ylab = "northing", zlab = "Elevation", clab = "elevation (m a.s.l.)",
	#        expand = 0.09, d = 2, phi = 25, theta = 40, resfac = 5,
	#        colkey = list(side = 2, length = 0.25))
	#dev.off()
}