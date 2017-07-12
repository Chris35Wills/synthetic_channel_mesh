# Program: 14c_get_minimum_surface_from_stack.r
# Create a new raster just of the minimum values of all single channel rasters, compiled into a stack
# Min surface should elimnate ridges where channels overlap
#
# Help and reading:
# http://www.inside-r.org/node/93677
# http://www.inside-r.org/packages/cran/raster/docs/stackSelect
# https://stat.ethz.ch/pipermail/r-sig-geo/2013-October/019572.html
#
#@author Chris Williams
#@date 27/04/16

library(raster)

check_rasters_for_all_zeros<-function(ras_path, ras_glob_name){
	
	rasList=Sys.glob(paste0(ras_path, ras_glob_name))
	empty_rasters=data.frame(paths_to_ignore=character())

	for (ras in rasList){
		
		a_raster=raster(ras)
		if (minValue(a_raster) == maxValue(a_raster)){  ### change this to a range check... omit NA
			empty_rasters=rbind(empty_rasters, ras) # append to list
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

get_min_from_stack<-function(raster_path_list, opath){
	
	ras_stack=stack(raster_path_list) #create raster stack

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

#################
## IMPLEMENTATION

### 79N 
#ras_path="C:/analysis_outputs_TEMP/Godthabsfjord_SYNTH_output/out_I_and_II_second_attempt/"
#ras_glob_name="channel_point2ras_agregated_channel_*_res_*m.tif"
#opath=ras_path

#get_min_from_stack(ras_path, ras_glob_name, opath)

## OLDER

### Godthabsfjord
ras_path="C:/analysis_outputs_TEMP/Godthabsfjord_SYNTH_output/out_I_and_II/"
ras_glob_name="channel_point2ras_agregated_channel_*.tif"
opath="C:/analysis_outputs_TEMP/Godthabsfjord_SYNTH_output/out_I_and_II/"

get_min_from_stack(ras_path, ras_glob_name, opath)


#set paths
	#tif_path="C:/GitHub/synthetic_channels/test_output/site_6_meshing_sep_channels/"
	#opath="C:/GitHub/synthetic_channels/test_output/site_6_meshing_sep_channels/"
#open rasters
	#tif=c((capture.output(cat(tif_path, "test_spline_150_m_1200m_every_20_channel_1.tif", sep=""))),
	#	  (capture.output(cat(tif_path, "test_spline_150_m_1200m_every_20_channel_2.tif", sep=""))),
	#	  (capture.output(cat(tif_path, "test_spline_150_m_1200m_every_20_channel_3.tif", sep=""))),
	#	  (capture.output(cat(tif_path, "test_spline_150_m_1200m_every_20_channel_4.tif", sep=""))),
	#	  (capture.output(cat(tif_path, "test_spline_150_m_1200m_every_20_channel_5.tif", sep=""))))

#### compare to normal raster
#
#min_path="C:/GitHub/synthetic_channels/test_output/site_6_meshing_sep_channels/"
#min=c((capture.output(cat(min_path, "test_spline_150_m_1200m_every_20_.tif", sep=""))))
#min_dem<-raster(min[1])
#
#old_path="C:/GitHub/synthetic_channels/test_output/site_6_meshing/"
#old=c((capture.output(cat(old_path, "test_spline_150_m_1200m_every_20_.tif", sep=""))))
#old<-raster(old[1])
#
#coords<-"+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#projection(min_dem) <- coords
#projection(old) <- coords
#
#uc <- overlay(old, min_path, fun=function(x,y){return(x+y)})
#
#fun <- function(r1, r2) { r1 - r2 }
#rc1 <- calc(old, min_dem, fun)
#
#diff<-old-min_dem # positive values where old dem was higher... expected at the inter-channel ridges
#image(diff)
#
#spline_xyz<-rasterToPoints(spline_rst)
#colnames(spline_xyz)<-c('x','y','z')
#spline_xyz_df<-as.data.frame(spline_xyz)
#
###raster z values as matrix (rather than raster)
#zData <-round(as.matrix(spline_rst),1)
#
#	
##ofile_name=capture.output(cat("/scratch/glacio1/cw14910/synthetic_channels_SCRATCH/better_dem_3d/test_spline_250_m_1200m_3D_", angle_string, "_deg.png", sep=""))
#png(ofile_name)
#angle=120
#persp3D(z=zData, xlab = "longitude", bty = "g",
#  	  ylab = "latitude", zlab = "Elevation", clab = "elevation (m a.s.l.)",
#        expand = 0.09, d = 2, phi = 15, theta = angle, resfac = 5,
#        colkey = list(side = 2, length = 0.25))
#dev.off()
