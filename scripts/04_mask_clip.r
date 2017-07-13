library(raster)
library(sp)
library(ggplot2)
library(tools)

# Plot points over raster
plot_mask_and_points<-function(mask, mask_extent=extent(-10000, 25000, -750000, -700000), pnts, centre=TRUE){
	
	mask_crop=crop(mask, mask_extent)
	map.p <- rasterToPoints(mask_crop)
    map.df <- data.frame(map.p)
    colnames(map.df)<-c('x','y','Mask')

	if (centre==TRUE){
		mask_pnt_plot=ggplot(map.df, aes(x,y)) +
      		geom_raster(aes(fill=Mask)) +
      		geom_point(data=pnts, aes(cl_x,cl_y, colour='red'), pch=20, na.rm=TRUE)
    } else {
	    mask_pnt_plot=ggplot(map.df, aes(x,y)) +
      		geom_raster(aes(fill=Mask)) +
      		geom_point(data=pnts, aes(x,y, colour='red'), pch=20, na.rm=TRUE)
    }

    mask_pnt_plot = mask_pnt_plot + 
    				theme_bw() +
      				coord_equal() +
      				xlab("Easting") +
      				ylab("Northing")

    mask_pnt_plot 
}

# The normal points go eitehr side of the centreline
# The centreline is not always the true centre but will be someoewhere near
# Consequently, the lengths of the normals either side should be approx equal
# This is not always the case such as where a channel passes a confluence - in 
# 	this case the normal can go up the confluence before it reaches non-ocean 
# 	portions of the mask which clip it
#
# Where a normal on one side is more than 1.5x the length of the other, it is 
# 	cropped to the max length of its opposite 
#   e.g. right side = 1000 m
# 		 left side  = 2000 m
# 		 corrected left side = 1000 m
clip_side_length_outliers<-function(data, diff_factor=1.5, verbose=0){

	length_before=nrow(data)
	left=subset(data, side==2)
	right=subset(data, side==1)

	# If dist of node from centre is greater than the diff_factor multiplied by the max dist of the other side, 
	#	limit it to the max length of the other side multiplied by the diff_factor
	if (max(right$cl_dist_r) > (diff_factor*max(left$cl_dist_l))) {
		right=right[right$cl_dist_r<=(diff_factor*max(left$cl_dist_l)),]
	} else if (max(left$cl_dist_l) > (diff_factor*max(right$cl_dist_r))) {
		left=left[left$cl_dist_l<=(diff_factor*max(right$cl_dist_r)),]
	}	

	data=rbind(left, right)
	
	length_after=nrow(data)

	if (verbose!=0){
		print(paste0("Length before: ", length_before))
		print(paste0("Length after: ", length_after))
	}

	return(data)

}

#' Keep only normals that are less than or equal to the mean normal length
#' we add the mean to the mean to give some more leeqay i.e. to relax the 
#' clipping a bit - otherwise this ovely favours the shorter segments....
second_pass_normal_filter<-function(data){
	
	print("Running second filter pass....")
	print(paste0("Length before: ", nrow(data)))
	#data=all_keeps

	m_length_l=mean(data$cl_dist_l)
	m_length_r=mean(data$cl_dist_r)

	limit_length=min(m_length_l, m_length_r)

	data=data[data$cl_dist_l<=(limit_length+(limit_length)),] # keep only normals that are less than or equal to the mean normal length
	data=data[data$cl_dist_r<=(limit_length+(limit_length)),] # we add the mean to the mean to give some more leeqay i.e. to relax the clipping a bit - otherwise this ovely favours the shorter segments....

	print(paste0("Length after: ", nrow(data)))
	
	return(data)

}

#' Clip mesh according to mask with OPTIONAL additional clip factor to minimize overall normal length 
#' This optional length clip comes into play where a channel goes past a confluence and there is potential 
#' for the normal to go up another channels centreline
mask_clipper<-function(path, mask, glob_path="*REARRANGED.csv", limit_different_edge_lengths=FALSE, diff_factor=1.25){
	
	filenames=Sys.glob(paste0(path, glob_path))

	if(length(filenames)==0){
		stop("No files found. Check you glob string (glob_path).")
	}

	for (f in filenames){
		#f=filenames[1]

		print(paste0("Working on: ", basename(f)))

		data=read.csv(f)
		coordinates(data)<-~x+y
		data$mask_value=extract(mask, data)
		data=as.data.frame(data)

		xmin=min(data$x)-4000
		xmax=max(data$x)+4000
		ymin=min(data$y)-4000
		ymax=max(data$y)+4000

		plotting=FALSE

		if (plotting==TRUE){plot_mask_and_points(mask, mask_extent=extent(xmin, xmax, ymin, ymax), pnts=data, centre=FALSE)}

		# Drop points if not ocean
		dropped_data=data[(data$mask_value!=0),]
		data=data[!(data$mask_value!=0),] # drop all points where mask is not 0

		if (plotting==TRUE){plot_mask_and_points(mask, mask_extent=extent(xmin, xmax, ymin, ymax), pnts=data, centre=FALSE)}

		# Drop points where there is a break in the cumualtive distance from centreline
		# 	--	delete rows from remaining data frame that have larger distances from the centreline than points already removed sharing the same centreline point

		all_fids=unique(data$cl_id)
		all_keeps=data.frame('x'=NA,
							 'y'=NA,
							 'cl_x'=NA,
							 'cl_y'=NA,
							 'cl_id'=NA,
							 'cl_dist_r'=NA,
							 'cl_dist_l'=NA,
							 'path'=NA,
							 'side'=NA,
							 'mask_value'=NA)

		for (fid in all_fids){
			
			#print(paste0("Working on: ", fid))
			
			# subset keeps by fid and select left and right 
			keeps=subset(data, cl_id==fid)
			keeps_r=subset(keeps, side==1)
			keeps_l=subset(keeps, side==2)

			# subset drops by fid and select left and right 
			drops=subset(dropped_data, cl_id==fid)
			drops_r=subset(drops, side==1)
			drops_l=subset(drops, side==2)

			# calculate the minimum cl_dist for this fid that was dropped 
			min_dropped_cl_dist_r=min(drops_r$cl_dist_r)
			min_dropped_cl_dist_l=min(drops_l$cl_dist_l)

			# Keep only remianing points less than the min cl_dist value (from the drops) 
			# Anything greater in the keeps dataframe is in another channel and we need to get rid of it

		#	print(paste0("Points before second drop (right): ", nrow(keeps_r)))
		#	print(paste0("Points before second drop (left): ", nrow(keeps_l)))
			keeps_r=keeps_r[keeps_r$cl_dist_r<min_dropped_cl_dist_r,]
			keeps_l=keeps_l[keeps_l$cl_dist_l<min_dropped_cl_dist_l,]
		#	print(paste0("Points after second drop (right): ", nrow(keeps_r)))
		#	print(paste0("Points after second drop (left): ", nrow(keeps_l)))

			keeps=rbind(keeps_r, keeps_l)
			#x11()
			#plot_mask_and_points(mask, mask_extent=extent(min(keeps$x)-5000, max(keeps$x)+5000, min(keeps$y)-5000, max(keeps$y)+5000), pnts=keeps, centre=FALSE)


			if (limit_different_edge_lengths==TRUE){
				keeps=clip_side_length_outliers(keeps, diff_factor=diff_factor, verbose=0)
			}

			#x11()
			#plot_mask_and_points(mask, mask_extent=extent(min(keeps$x)-5000, max(keeps$x)+5000, min(keeps$y)-5000, max(keeps$y)+5000), pnts=keeps, centre=FALSE)

			# append r and l to a dataframe....
			all_keeps=rbind(all_keeps,keeps)
		}

		all_keeps=na.omit(all_keeps)
	
		if (limit_different_edge_lengths==TRUE){
			all_keeps=second_pass_normal_filter(all_keeps)
		}
		
		if (plotting==TRUE){
			plot_mask_and_points(mask, mask_extent=extent(xmin, xmax, ymin, ymax), pnts=all_keeps, centre=FALSE)
			}

		# Write to file
		ofile=paste0(file_path_sans_ext(f), "_CLIPPED.csv")
		write.csv(all_keeps, ofile, row.names=FALSE)
	}
		#x11()
		#print("PLOTTING SHOULD NOW TAKE PLACE...")
		#plot_mask_and_points(mask, mask_extent=extent(xmin, xmax, ymin, ymax), pnts=all_keeps, centre=FALSE)

	cat(" \n")
	cat("*********************\n")
	cat("Normal clip complete.\n")
	cat("*********************\n")
}


if (getOption('run.main', default=TRUE)) {
	print("Run from import ... now running code with example code (will fail if earlier scripts in the processing chain have not already been run)")

	# Test data
	path="../test_outputs/"
	maskF="../test_data/aoi_mask.tif"
	mask=raster(maskF)
	glob_path="*REARRANGED.csv"
	mask_clipper(path, mask, limit_different_edge_lengths=TRUE, diff_factor=2, glob_path=glob_path)
}