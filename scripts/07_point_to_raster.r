
# #########################################################
# Program: 08_spline_all_points.r
#
# Reads in all points for development of the complete DEM (fitting a single spline to everything*)
# 
# Only uses points of one channel at a time - each grid being saved separately
#
# * in future, this will more likely be modified in future merging an interpolation of the channels 
#   with the surface dem interpolation (each portion of the dem having been inyterpolated independently)
#
# @Chris Williams 01/04/16
#########################################################

library(ggplot2)
library(gstat)
library(sp)
library(maptools)
library(raster)
library(scales)
library(reshape2)
library(plot3D)
library(fields)
library(ascii)
library(rgdal)
library(tools)
library(akima)

#install.packages("ggplot2")
#install.packages("gstat")
#install.packages("sp")
#install.packages("maptools")
#install.packages("raster")
#install.packages("scales")
#install.packages("reshape2")
#install.packages("plot3D")
#install.packages("rgeos")

##################
## FUNCTIONS
##################

#' Subsamples to every nth (n) element of a dataframe (dataframe)
Nth.row <- function(dataframe, n) dataframe[(seq (n, to=nrow(dataframe), by=n)),]

#' Takes in gridded points, aggregates them (lower resoluton) by a factor (aggFactor) and saves them as a raster (ras_out)
#' Default aggregation function is mean
aggregate_and_save<-function(point2ras, aggFactor, ras_out){
	# Aggregate 
	point2rasAGG=aggregate(point2ras, fact=aggFactor, fun=mean, na.rm=TRUE)
	# Save as raster
	cat("\nSaving raster...\n")
	writeRaster(point2rasAGG, ras_out, overwrite=TRUE)
	cat("\n Raster saved:", ras_out)
}

##################
## TESTS
##################

#' Tests to see if all column names are equal before proceeding
test_col_names<-function(list_1, list_2, list_3){
	list_1_names=colnames(list_1)	
	list_2_names=colnames(list_2)
	list_3_names=colnames(list_3)

	if(identical(list_1_names, list_2_names)==FALSE){
	stop("\nHeaders of all xyz versions of input files must match as x y z\n")
	}
	else if(identical(list_1_names, list_3_names)==FALSE){
	stop("\nHeaders of all xyz versions of input files must match as x y z\n")
	}
	else{
		# do nothing
	}
}

##################
##########

##################
## MAIN CODE
##################

#' @param aggFactors = aggregation factor - value of 1 returns channel raster at resolution of the mask
points_to_raster<-function(path, mask_path, mask, opath, glob_extension="", aggFactors=c(1)){
	
	
	cat("Opening data...\n")
	
	if (glob_extension!=""){
		print("Finding channels...")
		channels=Sys.glob(paste0(path, glob_extension))
	} else if (glob_extension==""){
		stop("No channels provided and not sure what extension to use to find them.")
	}

	maskF=c(capture.output(cat(mask_path, mask, sep="")))
	mask=raster(maskF)

	##################
	##Extract xyz

	cat("\nExtracting xyz...\n")
	for(i in 1:length(channels)){
	#	i=1

		this_channel=channels[i]
		split1=strsplit(this_channel, "path_")[[1]][2]
		channel_num=strsplit(split1, "_clipped")[[1]][1]
		
		#create empty raster same dimensions as land dem
		#https://cran.r-project.org/web/packages/raster/vignettes/Raster.pdf
		emptyGrid=raster()
		res(emptyGrid)[1]<-res(mask)[1]
		res(emptyGrid)[2]<-res(mask)[2]
		extent(emptyGrid)<-mask@extent
		crs(emptyGrid)<-mask@crs
		nrow(emptyGrid)<-mask@nrows
		ncol(emptyGrid)<-mask@ncols

		compareRaster(emptyGrid, mask)

		# getting xyz data
		cat("\nExtracting normal xyz...\n")
		channel=na.omit(as.data.frame(read.csv(channels[i])))

	#	channel_xy=channel[c('x','y')]
	#	channel_z=channel[,c('norm_elev')]
		channel_xyz=channel[,c('x','y','norm_elev')]
		colnames(channel_xyz)<-c('x','y','z')

		cat("\nExtracting centreline xyz...\n")
		centre_xyz=channel[,c('cl_x','cl_y','cl_elev')]
		centre_xyz=centre_xyz[!duplicated(centre_xyz), ]
		colnames(centre_xyz)<-c('x','y','z')
	#	centre_xy=centre_xyz[,c('x','y')]
	#	centre_z=centre_xyz[,c('z')]
		

		cat("\nCombining normal and centreline xyz...\n")
		xyz=rbind(channel_xyz, centre_xyz)
	#channel_xy=rbind(channel_xy, centre_xy)
	#channel_z=rbind(channel_z, centre_z)

		## Stop wasting memory :)
		channel_xyz=NULL

		##################
		## Create as dataframe
		xyz_df=data.frame(xyz)
		xyz=NULL

	#point2ras=rasterize(channel_xy, emptyGrid, channel_z, fun=mean)
		point2ras=rasterize(xyz_df[c('x','y')], emptyGrid, xyz_df$z, fun=mean)

		##################
		## Aggregate to lower resolution
		
		#aggFactors=c(2,3,4)	#c(4)
		for (aggFactor in aggFactors){
			#ras_out= capture.output(cat(opath, "channel_point2ras_agregated_channel_",i,"_aggFactor_", aggFactor,".tif", sep=""))
			ras_out= capture.output(cat(opath, "channel_point2ras_agregated_channel_",channel_num,"_res_", aggFactor*res(mask)[1], "m.tif", sep=""))

			if (aggFactor>1){
				aggregate_and_save(point2ras, aggFactor, ras_out)	
			} else {
				cat("\nSaving raster...\n")
				writeRaster(point2ras, ras_out, overwrite=TRUE)
				cat("\nRaster saved:", ras_out)
			}
		}

		image(point2ras)
			
		xyz_df=NULL
		
	}

	cat("Check coordinate systems of each gdalinfo *.tif\n")
	cat("Convert otherwise using: gdal_translate -a_srs '+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' JUST_CHANNE_spline_150_m_1200m_every_1_channel_1.tif JUST_CHANNE_spline_150_m_1200m_every_1_channel_1_PROJ.tif\n")
	cat("Created and saved all surfaces - END\n")

}

#' Resample raster and interpolate using a bicubic spline
#' 
#' raster_paths  - a list of raster paths - routin will be run on each raster
#' new_cell_size - cell size to resample to
#' maskPath      - a raster the extent of which you wish the output to match
resample_rasters<-function(raster_paths, new_cell_size, maskPath){
	
	for (ras1 in raster_paths){

		#https://rdrr.io/cran/akima/man/interp.html

		print(paste0("Resampling and bilinear interpolation of", basename(ras1), " to ", round(new_cell_size), "m..."))


		outF_resampled=paste0(file_path_sans_ext(ras1), "_RESAMPLED_", round(new_cell_size), "m.tif")
		outF_bc=paste0(file_path_sans_ext(ras1), "_RESAMPLED_", round(new_cell_size), "m_BC.tif")

		ras1=raster(ras1)
		xyz=as.data.frame(rasterToPoints(ras1))
		colnames(xyz)=c('x','y','z')

		# Create target grid - extent of points for now (will resamp to main AOI below)
		new_extent=extent(min(xyz$x),max(xyz$x),min(xyz$y),max(xyz$y))
		targetGrid=raster()
		res(targetGrid)[1]<-new_cell_size
		res(targetGrid)[2]<-new_cell_size
		extent(targetGrid)<-new_extent # ras1@extent 
		crs(targetGrid)<-ras1@crs
		nrow(targetGrid)<-round((new_extent@ymax-new_extent@ymin)/new_cell_size) # round((ras1@extent@ymax-ras1@extent@ymin)/new_cell_size)
		ncol(targetGrid)<-round((new_extent@xmax-new_extent@xmin)/new_cell_size) # round((ras1@extent@xmax-ras1@extent@xmin)/new_cell_size)


		xy_target=as.data.frame(coordinates(targetGrid))

		### Now Akima splines with accurracy of bicubic polynomial for irregular gridded data
		print("Applying bicubic polynomial...")
		iRspl <- with(xyz,interp(x,y,z,linear=FALSE,nx=nrow(targetGrid),ny=ncol(targetGrid)))
			## PLOT
			# breaks <- pretty(c(min(xyz$z),max(xyz$z)),10)
			# colors <- heat.colors(length(breaks)-1)
			# image(iRspl,breaks=breaks,col = colors)
			# contour(iRspl,col="black",levels=breaks,add=TRUE)

		# convert bc to raster
		bc_raster=raster(iRspl)
		crs(bc_raster)<-ras1@crs

		# omit bc where original raster was na
		bilinear_resample=resample(ras1, targetGrid, method="bilinear")
		na_mask=!is.na(bilinear_resample) # 1 if not NA - without the '!' this would be the other way around
		na_mask[na_mask==0]=NA
		bc_raster_resamp=resample(bc_raster, bilinear_resample)
		compareRaster(bc_raster_resamp, bilinear_resample)
		bc_raster_masked=bc_raster_resamp*na_mask

			## PLOT
			# image(bc_raster_masked,breaks=breaks,col = colors)
			# contour(bc_raster_masked,col="black",levels=breaks,add=TRUE)

		#resamp 2 (resamp to AOI grid dimensions - required later for stacking)
		#print("Resampling to AOI grid dimensions...")
		mask=raster(maskPath)
		resamp_raster_masked<-resample(bc_raster_masked, mask, method="ngb")
		
		# write to file
		#writeRaster(bc_raster_masked, outF_bc, overwrite=TRUE)
		writeRaster(resamp_raster_masked, outF_bc, overwrite=TRUE)


	}

}

################
# IMPLEMENTATION

## Test AOI
# 
mask_path="../test_data/"
path="../test_outputs/"
glob_extension="*___PIECEWISE.csv"
opath=path
aoi_mask=paste0(mask_path, "godthabsfjord_mask__CROP.tif")

points_to_raster(path, mask_path, aoi_mask, opath, glob_extension=glob_extension, aggFactors=c(1))


## resample_rasters(raster_paths=Sys.glob(paste0(path, glob_extension), 
##					new_cell_size=500, 
##					maskPath=aoi_mask)