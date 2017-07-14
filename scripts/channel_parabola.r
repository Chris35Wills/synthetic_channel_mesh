#source("C:/GitHub/R_bits/my_functions/plotting.r")
#source("C:/Github/synthetic_channels/meshing/sorted_scripts/parabola_funcs_where_obs_available.r")        
#source("C:/GitHub/synthetic_channels/meshing/sorted_scripts/channel_parabola.r")

#library(fields) # for splies (er.g. Tps)
#library(ggplot2)
#library(rgl) # for 3d plotting
library(raster)
library(sp)
#library(car)
#source("./plotting.r")

source("./parabola_funcs_where_obs_available.r")        


read_me_df<-function(data_f){
  data=read.csv(data_f)
  df<-as.data.frame(data)
  return(df)
}



# Provides test data
# as_df  - is TRUE returns a list containing dataframes, otherwise just a list of file path strings
# sparse - 2 - data available along length of fjord
#        - 1 - data available along parts of fjord
#        - 0 - no data available in fjord (relies on the seed and mouth input points)
open_test_data<-function(data_path="C:/GitHub/synthetic_channels/meshing/test_data/fenty_aoi_region/", as_df=TRUE, sparse=2){

  print("Opening test data")
  normal_points_F=capture.output(cat(data_path, "test__densified_path_1_clipped_normals_avg_6_pnts_1200m_REARRANGED_CLIPPED___fenty_aoi.csv", sep=""))      # FORMAT: AS PER *.csv
  edge_elevations_F=capture.output(cat(data_path, "test_densified_path_1_clipped_normals_avg_6_pnts_1200m_REARRANGED_CLIPPED_edge_elevs___fenty_aoi.csv", sep="")) # FORMAT: AS PER *.csv
  
  if (sparse==2){
    observations_F=capture.output(cat(data_path, "test_OMG_obs_centreline_xyz__fenty_aoi.csv", sep=""))  # FORMAT: X,Y,Z
  } else if (sparse==1){
    observations_F=capture.output(cat(data_path, "test_OMG_obs_centreline_xyz__fenty_aoi__SPARSE.csv", sep=""))  # FORMAT: X,Y,Z
  } else if (sparse==0){
    observations_F=capture.output(cat(data_path, "test_OMG_obs_centreline_xyz__fenty_aoi__NO_OBSERVATIONS.csv", sep=""))  # FORMAT: X,Y,Z
  }
  
  f_out=capture.output(cat(data_path, "test_spline_fjord___with_observations__PIECEWISE_LINEAR.csv", sep=""))

  if (as_df==TRUE){
    df_norm=read_me_df(normal_points_F)
    df_edge=read_me_df(edge_elevations_F)
    df_obs=read_me_df(observations_F)
  
    return(list(df_norm, df_edge, df_obs, f_out))  
  } else {
  
    return(list(normal_points_F, edge_elevations_F, observations_F, f_out))  
  }

}

#' Calculate unknowns in parabola equation
calc_parabola_vertex <- function(x1, y1, x2, y2, x3, y3) {
  #Pinched and subsequently modifed to get the unknowns for defining the parabola:
  #http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
  denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
  A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
  B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom
  C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom
  abc<-c(A,B,C)
  return(abc)
}


#' Uses the 'fields' library
#' Create spline using x, y and z - returns a surface
spline_surface<-function(xyz, plotting=T){

  combo=as.data.frame(xyz)

  fit<- Tps(combo$x, combo$y) 

  coord<-data.frame(combo$x, combo$y)
  z<-data.frame(combo$z)

  fit<-Tps(coord,z)

  xg<-make.surface.grid(fields.x.to.grid(coord)) # makes a mesh grid of x and y
  if (plotting==TRUE){plot(xg)} # displays the mesh grid
  fhat<- predict( fit, xg) # indexes values of based on the thin plate spline to their xy index locations (i.e. this is 1d and is therefore a length of the dimensions of plot(xg))
  out.p<- as.surface( xg, fhat)
  
  return(out.p)
}


#' Tests if 2 variables are equal to one another - exits if not
#'
#' @param a some input (int, double, float, string...) 
#' @param b some input (int, double, float, string...)
#' @examples
#' expect_equal("test", "test1")
#' Error in expect_equal("test", "test1") :  1st input is not equal to 2nd input
expect_equal<-function(a,b){
  if(a!=b){
    stop("1st input is not equal to 2nd input")
    warning("1st input is not equal to 2nd input")
  }
}


#' Tests if a variable is equal to 0 - exits if not
#'
#' @param a some input (int, double, float, string...) 
#' @examples
#' val<-1
#' expect_zero(val)
#' Error in expect_zero(val) : Input is not equal to zero
expect_zero<-function(a){
  if(a!=0){
    stop("Input is not equal to zero")
    warning("Input is not equal to zero")
  }
}


expect_matching_fid<-function(mx_1, mn_1, mx_2, mn_2){
  
  if(mx_1!=mx_2){
    stop("Max FIDs don't match between the input normals and the edge points")
    warning("Check you have the right normal/edge elevation file pair")
  }
  
  else if(mn_1!=mn_2){
    stop("Min FIDs don't match between the input normals and the edge points")
    warning("Check you have the right normal/edge elevation file pair")
  }
  
}


expect_negative_slope<-function(seed_elevation,mouth_elevation){
  if(as.numeric(seed_elevation) < as.numeric(mouth_elevation)){
    stop("Seed is lower than the mouth elevation - this will cause an upward sloping profile")
    warning("To allow for an upward sloping profile, you'll need to disable this catch.")
  }
}

# centre_dist - if edges are equal distance from the centre, this will be zero
test_parabola_func<-function(a,b,c,seed_elevation, centre_dist=0){
    #centre_dist=0 # assumes the centre node is truely between the 2 edges (which are at positive and negative distances e.g. -2000 and 2000)
    r_depth = (a*(centre_dist^2))+(b*centre_dist)+c
   
    # Base comparison on values rounded to nearest 0.1 (otherwise get rounding errors...)
    r_depth=round(r_depth*10)/10
    seed_elevation=round(seed_elevation*10)/10

    #print(paste("r_depth: ", r_depth))
    #print(paste("seed_elevation: ", seed_elevation))

    if (is.na(r_depth) | is.na(seed_elevation)){
      cat("*******************\n")
      cat("*******************\n")
      print("Testing parabola func")
      print(paste0("Your seed elevation was: ", seed_elevation))
      print(paste0("Your depth was: ", r_depth))
      #stop("Testing parabola func\nYou have na values...\nIf using a model to fit cl elevations, check you are in the same range as the model")
      cat("Testing parabola func\nYou have na values...\nIf using a model to fit cl elevations, check you are in the same range as the model")
      cat("*******************\n")
      cat("*******************\n")
      return(parabola_seed_depth_non_match=TRUE)
    } else if(isTRUE(all.equal(r_depth, seed_elevation))==FALSE){ # safest way to see if small numbers are equal (better than using a == b)
      cat("*******************\n")
      cat("*******************\n")
      cat("Calculated elevation:", r_depth, "\n")
      cat("Centre elevation:", seed_elevation, "\n")
      #stop("Centre of parabola not equal to the seed elevation (which is used as the centreline elevation - check call to calc_parabola_vertex()")
      cat("Centre of parabola not equal to the seed elevation (which is used as the centreline elevation - check call to calc_parabola_vertex()")
      cat("*******************\n")
      cat("*******************\n")
      return(parabola_seed_depth_non_match=TRUE)
    } else {
      return(parabola_seed_depth_non_match=FALSE)
    }

}

test_make_parabola<-function(){
  
  dist_from_centre=seq(from=-400, to=400, by=10)

  edge_1=c(-400,20) # (dist from centre, elevation)
  centre=c(0,-200)  # (dist from centre, elevation)
  edge_2=c(400,100) # (dist from centre, elevation)

  abc<-calc_parabola_vertex(edge_1[1], edge_1[2], centre[1], centre[2], edge_2[1], edge_2[2])
  a=abc[1]
  b=abc[2]
  c=abc[3] # mid point

  ##################
  #Define points along parabola (according to distance from centre)

  row_est=length(dist_from_centre)
  norm_elevs <- matrix(NA, nrow=row_est, ncol=1) # pre-define list size to hold values
  norm_cl_dists <- matrix(NA, nrow=row_est, ncol=1) # pre-define list size to hold values

  #norm_r_elev
  for (x in 1:length(dist_from_centre)){ 
    
    x_val=dist_from_centre[x] # <<< centreline distance
    
    elev = (a*(x_val^2))+(b*x_val)+c
    #elev = round(elev)

    #cat("norm elev: ", elev, "\n")

    norm_elevs[x]<-rbind(elev)
    norm_cl_dists[x]<-rbind(x_val)
  }

  norm_elevs=na.omit(norm_elevs)[,1]
  norm_cl_dists=na.omit(norm_cl_dists)[,1]

  cat("Known point:\n")
  cat("Point 1 : -400, 20\n")
  cat("Point 2 : 0, 20\n")
  cat("Point 3 : 400, 100\n")

  norm_df=data.frame("norm_cl_dists"=norm_cl_dists, "norm_elevs" = norm_elevs)
  ggplot(norm_df, aes(norm_cl_dists, norm_elevs), colour='black') +
    geom_point() +
    geom_point(data=subset(norm_df, norm_cl_dists == edge_1[1]), aes(norm_cl_dists, norm_elevs), size=3, colour='green') +
    geom_point(data=subset(norm_df, norm_cl_dists == edge_2[1]), aes(norm_cl_dists, norm_elevs), size=3, colour='red') +
    geom_point(data=subset(norm_df, norm_cl_dists == centre[1]), aes(norm_cl_dists, norm_elevs), size=3, colour='orange')

  #plot(norm_cl_dists,norm_elevs)
}


plot_test_output<-function(norms_csv){
  
    print("Attempting to plot test output...")

    norm_data=read.csv(norms_csv)
    df<-as.data.frame(norm_data)
  
    # 2d plot
    p1<-ggplot(df, aes(x, y)) + 
      #geom_point(aes(colour=norm_elev)) + 
      geom_point()+
      geom_point(data = df, aes(cl_x, cl_y, colour=cl_elev)) +
      coord_fixed() +
      ggtitle("Normal and centreline locations...")

    p2<-ggplot(df, aes(cl_id, norm_elev)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      ggtitle("Linear model (not necessarily\nused to calc elevs)")

    p3<-ggplot(df, aes(cl_id, cl_elev)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      ggtitle("Linear model (not necessarily\nused to calc elevs)")
    
    multiplot(p1, p2, p3, cols=3)

}


#' true_centre : if TRUE, then the centreline elevations (when defining the parabola) are assigned to the "true centre" i.e. which 
#'               lies exactly halfway between the edges of a given normal. Otherwise the centreline elevations are assigned to the 
#'               mapped centreline position. Using one or the other will result in slightly different normal node elevations as you 
#'               shift the position of one of your "known" values used to define the parabola.
#'
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1)  # piecewise
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1) # piecewise
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0)  # linear
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0) # linear
#
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1, sparse=1, dist_threshold=2000)  # piecewise + sparse
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1, sparse=1, dist_threshold=2000) # piecewise + sparse
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0, sparse=1)  # linear + sparse
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0, sparse=1) # linear + sparse
#
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1, sparse=0)  # piecewise + no obs
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1, sparse=0) # piecewise + no obs  - OK
#test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0, sparse=0)  # linear + no obs     - OK
#test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0, sparse=0) # linear + no obs     - OK
#
# WHAT TO EXPECT from the tests
#
# Where true_centre==TRUE, the observations will be slightly offset from expected values - this is because the expected values are 
#   calculated from a model using th exact centre of the channel - elevations are then calculated using the nodes along each normal 
#   which are in a sightly different position which alters the elevations calculated for them and the caccumulated distance between 
#   them
# Where true_centre==FALSE, observations should match the expected as the same nodes and therefore the same accumulated distance are 
# used to inform the model predicting centreline elevations
# Setting the dist_threshold too small will result in observtions not being assigned to normals along the centreline
test_channel_parabola__piecewise<-function(data_path="C:/GitHub/synthetic_channels/meshing/test_data/fenty_aoi_region/", true_centre=TRUE, piecewise=1, sparse=2, dist_threshold=2000){
  
  dat=open_test_data(data_path, as_df=FALSE, sparse=sparse)
  norm_f=dat[[1]]
  edge_f=dat[[2]]
  obs_f=dat[[3]]
  f_out=dat[[4]]
  
  seed_xyz=data.frame(x=-717233,y=-1326743,z=-158)
  mouth_xyz=data.frame(x=-700582,y=-1358375,z=-920)

  dat_out=channel_parabola__piecewise(norm_f, edge_f, obs_f, f_out, seed_xyz, mouth_xyz, true_centre=true_centre, piecewise=piecewise, dist_threshold=dist_threshold)
  
  elev_mesh=as.data.frame(dat_out[1]) # normal and centreline points with assigned elevations
  xyz=as.data.frame(dat_out[2]) # points used to define the piece-wise model - if true_centre==TRUE, then these will be the EXACT centre (x,y) (not the nearest mesh node equvalents), otherwise they will be the mapped cl(x,y) positions


  ######################################
  ## EVERYTHING BELOW IS FOR PLOTTING ##
  ######################################

  if (true_centre==TRUE){
    if (piecewise==1){   
      if (sparse==2){
        opath=paste0(data_path, "plots_OBS_true_centre__piecewise/")
      } else if (sparse==1){
        opath=paste0(data_path, "plots_OBS_true_centre__piecewise__SPARSE/")
      } else if (sparse==0){
        opath=paste0(data_path, "plots_OBS_true_centre__piecewise__NO_OBS/")
      }
    } else {
      if (sparse==2){
        opath=paste0(data_path, "plots_OBS_true_centre__linear/")
      } else if (sparse==1){
        opath=paste0(data_path, "plots_OBS_true_centre__linear__SPARSE/")
      } else if (sparse==0){
        opath=paste0(data_path, "plots_OBS_true_centre__linear__NO_OBS/")
      }
    }   
  } else if (true_centre==FALSE){  
    if (piecewise==1){
      if (sparse==2){
        opath=paste0(data_path, "plots_OBS_mapped_centre__piecewise/")
      } else if (sparse==1){
        opath=paste0(data_path, "plots_OBS_mapped_centre__piecewise__SPARSE/")
      } else if (sparse==0){
        opath=paste0(data_path, "plots_OBS_mapped_centre__piecewise__NO_OBS/")
      }
    } else {
      if (sparse==2){
        opath=paste0(data_path, "plots_OBS_mapped_centre__linear/")
      } else if (sparse==1){
        opath=paste0(data_path, "plots_OBS_mapped_centre__linear__SPARSE/")
      } else if (sparse==0){
        opath=paste0(data_path, "plots_OBS_mapped_centre__linear__NO_OBS/")
      }
    }
  }

  dir.create(opath, recursive=TRUE)

  #plot the parabola for a given fid
  for (fid in min(elev_mesh$cl_id):max(elev_mesh$cl_id)){
     
    elevs_sample=subset(elev_mesh, elev_mesh$cl_id==fid)
    elevs_sample_LEFT=subset(elevs_sample, elevs_sample$side==2)
    elevs_sample_RIGHT=subset(elevs_sample, elevs_sample$side==1)
    centre_point=data.frame(dist=0, cl_elev=unique(elevs_sample_RIGHT$cl_elev)) # NOT NECESSARILY IN THE MIDDLE - has an x value of 0 on the axis for which the parabola was calculated

    elevs_sample_LEFT$type="Left normals"
    elevs_sample_RIGHT$type="Right normals"
    centre_point$type="Mapped Centre (0)"

    # create output plot file
    FORMATC <- function(x) formatC(x, width = 3,flag = 0) # pad numbers
    ofile=paste0(opath, FORMATC(fid), "_parabola_plot_fid_", ".png")
     
    dir.create(dirname(ofile), recursive=TRUE)
    parab_plot=ggplot(elevs_sample_LEFT, aes(cl_dist_l, norm_elev, colour=type)) +
      geom_point() +
      geom_point(data=elevs_sample_RIGHT, aes(cl_dist_l*-1, norm_elev, colour=type)) +
      geom_point(data=centre_point, aes(dist, cl_elev, colour=type)) +
      ggtitle(paste0("FID: ", fid))
     
    ggsave(parab_plot, file=ofile)
     
   }

  # Plot piece-wise model fit
  ofile=paste0(opath, "linear_model_check.png")


  if (piecewise==1){
    piecewise_check_plot=check_obs_against_centreline_model(xyz, elev_mesh, seed_xyz, mouth_xyz, true_centre=true_centre, piecewise=piecewise, verbose=verbose)
    ggsave(piecewise_check_plot, file=ofile)
  } else { 
   # plot fit to linear model....
   linear_check_plot=check_obs_against_centreline_model(xyz, elev_mesh, seed_xyz, mouth_xyz, true_centre=true_centre, piecewise=piecewise, verbose=verbose)
   ggsave(linear_check_plot, file=ofile)
  }
    
}

#' Runs test_channel_parabola__piecewise() for various test cases
run_test_channel_parabola__piecewise<-function(){
  # Lots of obs
  test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1)  # piecewise
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1) # piecewise
  #test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0)  # linear
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0) # linear

  # Some obs
  test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1, sparse=1, dist_threshold=2000)  # piecewise + sparse
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1, sparse=1, dist_threshold=2000) # piecewise + sparse
  #test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0, sparse=1)  # linear + sparse
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0, sparse=1) # linear + sparse

  # No obs
  test_channel_parabola__piecewise(true_centre=TRUE, piecewise=1, sparse=0)  # piecewise + no obs
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=1, sparse=0) # piecewise + no obs
  #test_channel_parabola__piecewise(true_centre=TRUE, piecewise=0, sparse=0)  # linear + no obs
  #test_channel_parabola__piecewise(true_centre=FALSE, piecewise=0, sparse=0) # linear + no obs
}

#' Calculate elevations along a channel cross section based on knonw edge elevations and a known centreline elevation
#' Data expected in specific formats:
#' df: x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side
#' df_edge: x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side,bed_elev_m
#' count - pertains to the iteration (i.e. if 3 transects, count will go up to 3)
#' last_ind_counter - a list the length of the number of fids to be considered - keeps track of the last index of the output dataframe that was written to so to avoid overwriting
#' true_centre=TRUE/FALSE
#' 
#'      DETAILS
#'      true_centre value
#'      The "centreline" mapped by the skeletonisation process is not always exactly between the edge nodes for a given normal - as edge nodes are calculated spearately to the centreline. 
#'      Edge points are given a "distance from the centreline value" - as edges are not necessarily equidistanct from the centre node, these values can differ
#'      When calculating a parabola from 3 points, the edges and centre points x values do not need to change
#'      e.g.
#'      
#'      | point  | dist from centre | elev |
#'      ---
#'      | edge 1 | -400 | 20 |
#'      | mid    | 0    | -200 |
#'      | edge 2 | 200 | 50 |
#'      
#'      The "dist from centre" values in this instance do not matter - the x scale for the parabola is then -400 -> 200
#'      
#'      However, there is another method (channel_parabola__piecewise()) which uses edge points and the nearest obsevation in the fashion above to calculate a parabola 
#'      for a given cross-section. The middle of this parabola is then considered where the middle is calculated as being between edge 1 and edge 2, using the same x scale - for the example 
#'      above, the "true midpoint" is at 200-(200-(-400))/2
#'      From this parabola, the known points are the edges and the "true midpoint"
#'      Thus, the "true midpoint" does not have an x value of 0 but of -100. This is important when considering the parabola to fit for a given location! Thus when fitting the parabola with 
#'       these points the midpoint cannot be assumed to be zero and therefore must be corrected to represent its true x value. Only the midpoint should be corrected!
#'      
#'      2 cases:
#'      
#'        parabola from 3 points - no need to correct midpoints x ("midpoint" is an approximation and the scale is btween the edge point x values)
#'        parabola from edges and true centre - scale is between the edge point x values and the "true centre" lies equally between them
#'                                             - as "distance from mid point" is a field in the table of normal points, no row has a value of 0 (this is the midpoint!)
#'                                               therefore when adding the midpoint values to the final dataframe in set_elevation(), the centreline elevation is calculated 
#'                                               with an x value of 0 - this will not necessarily match the input centrepoint as the input does not account for the defined 
#'                                               spacing of points along the normal, so expect a small difference (if the offset of the true midpoint is cleanly divisible by 
#'                                               the spacing of points along the normal, then the values may well match)
#'                                               - to force the centreline elevs to match the input, the centreline elev is caclualted using the precise offset of the true centre point from 0
#'      
#' add_centre_offset_to_output=FALSE - if TRUE, the offset of the true centre from the mapped centre is included in the output dataframe
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' 
#' Func can be called like:
#' 
#'      count=0
#'      #' # make output
#'      N <- 1e4  # some magic number (possibly an overestimate) for pre-setting rows of dataframe
#'      out_df <- data.frame(x=rep(NA, N), y=rep(NA, N),cl_x=rep(NA, N),cl_y=rep(NA, N),cl_id=rep(NA, N), cl_dist_r=rep(NA, N), cl_dist_l=rep(NA, N), path=rep(NA, N), side=rep(NA, N),norm_elev=rep(NA, N),cl_elev=rep(NA, N),stringsAsFactors=FALSE) 
#'      last_indx_counter <- matrix(NA, nrow=N, ncol=1)
#'      for (fid in 1:10){
#'         out=set_elevations(df, df_edge, fid, out_df, count, last_indx_counter, max_depth=-200)  
#'         out_df=out$out_df
#'         count=out$count
#'         last_indx_counter=out$last_indx_counter
#'      }
set_elevations <- function(df, df_edge, fid, out_df, count, last_indx_counter, centreline_depth=-200, verbose=1, true_centre=FALSE, add_centre_offset_to_output=FALSE){
  
  if(verbose==1){ cat("Working on line FID: ", fid, "\n") }
    
    #subset df according to fid
    cl_eg=subset(df, cl_id == fid)
    if(verbose==1){
      cat("Length of subset: ", nrow(cl_eg), "\n")
    }
    
    if(nrow(cl_eg)<2){

      return(outputs=FALSE)
    
    } else {

    # keep only cols required 
    keeps <- c("x","y","cl_x","cl_y","cl_id","cl_dist_r","cl_dist_l","path","side")
    cl_eg=cl_eg[keeps]

    # get xy normal of normal points
    keeps <- c("x", "y")
    xy=cl_eg[keeps]


    
    ##################
    #Define known points
    right_max_cl_dist<-max(cl_eg[cl_eg$side==1,]$cl_dist_r)
    left_max_cl_dist<-max(cl_eg[cl_eg$side==2,]$cl_dist_l)
    
    ##****************************************************************************##
    ##**** SET BANK ELEVATION ACCORDING TO NAREST BED ELEVATION DEM POINTS... ****##
    ##****************************************************************************##
    
    #Read in the edge point file for this dataset - this will be the output of *_get_channel_bank_elevation
    #Get the elevation ( a new column added in 05b*.py) associated with the right/left max distance - use this  
    #value to set the bank elevation (currently this is hardwired as a value of 0)
    
    EDGE_cl_eg=subset(df_edge, cl_id == fid) # subsets the edge elev dataframe to the specific fid
    right_edge_elev=subset(EDGE_cl_eg$bed.elev.m., EDGE_cl_eg$side==1)
    left_edge_elev=subset(EDGE_cl_eg$bed.elev.m., EDGE_cl_eg$side==2)
    
    edge_r=c(right_max_cl_dist*-1,right_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #1
    edge_l=c(left_max_cl_dist,left_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #2
    
          # if(verbose==1){
          #   cat("Right edge elevation: ", right_edge_elev, "m\n")
          #   cat("Left edge elevation: ", left_edge_elev, "m\n")
          # }
       
    # ***********
    # *********************
    # ***************************************
    # The mapped centreline point is not always exactly between the two normals (as these can vary slightly in length according to cliping based on masking). 
    # The right and left normal points are at set distances of a mid point, whose distance is 0 m - if points normal to the centre are 
    # not equal e.g. right goes -2000 m and left goes 1000 m from the centrepoint, then the normal total length is 3000m and so the real 
    # centre point is at -500 - this is corrected here...
          
    #Calculate the true centre according to max normal lengths (this is essentially a catch should left and right be different lengths)
    #If normals either side are equal, this value will be zero
    dist_from_cl_of_mid<-edge_r[1]+(edge_l[1]-edge_r[1])/2

    if (true_centre==TRUE){
      mid_depth=c(dist_from_cl_of_mid,centreline_depth) # < known point #3 (centre) - distance component set dynamically should the left and right half widths be different
      if(verbose==1){
        cat("Distance of true mid-point from geometrically defiend cl: ", dist_from_cl_of_mid, "m \n")
        } 
    } else {
      mid_depth=c(0,centreline_depth)
    }

    
    ##################
    #Calculate unknowns in parabola equation

    # As we assume our centreline elevation actually occurs at the centre, we must pass the centre nodes distance as being between the left and right edges (see variable: dist_from_cl_of_mid)
    # The correction is added to shift the 3 known points relative to zero - the actual locations of these points does not change of course
    # The parabola will now be centred on 0 which is vital considering the intercepts calculated!!
    if (true_centre==TRUE){
      cat("TRUE centre implementation...\n")
      #correction=mid_depth[1]*-1
      #abc<-calc_parabola_vertex(mid_depth[1]+correction, mid_depth[2], edge_r[1]+correction, edge_r[2], edge_l[1]+correction, edge_l[2]) # centre dist set to zero + move edges by the true centre/centre offset so centre point is now in the middle
      #abc<-calc_parabola_vertex(mid_depth[1]+correction, mid_depth[2], edge_r[1], edge_r[2], edge_l[1], edge_l[2]) # centre dist set to zero + move edges by the true centre/centre offset so centre point is now in the middle
      abc<-calc_parabola_vertex(mid_depth[1], mid_depth[2], edge_r[1], edge_r[2], edge_l[1], edge_l[2]) # centre dist set to zero + move edges by the true centre/centre offset so centre point is now in the middle
    } else {
      abc<-calc_parabola_vertex(mid_depth[1], mid_depth[2], edge_r[1], edge_r[2], edge_l[1], edge_l[2]) # centre dist set to zero + move edges by the true centre/centre offset so centre point is now in the middle
    }

    a=abc[1]
    b=abc[2]
    c=abc[3] # mid point
 
    if (true_centre==TRUE){
      parabola_seed_depth_non_match=test_parabola_func(a,b,c,mid_depth[2], centre_dist=mid_depth[1])
    } else {
      parabola_seed_depth_non_match=test_parabola_func(a,b,c,mid_depth[2]) # will fail if the mid_depth is not calculated at 0
    }
       
    if (parabola_seed_depth_non_match==FALSE){ # i.e. if test_parabola_func() passes
      ##################
      # Define points along parabola (according to distance from centre)
      # 
      # Each normal point has an associated distance from the centreline (column cl_dist_l or cl_dist_r - both the same..)
      # This distance is used as the x component of the parabola
      # Both left and right normals of the centre should be equal however as this is not the case, it must be accounted 
      # for using the correction factor calculated above    
      # When referring to the end output dataframe, by solely using the distances of points either side of the centre, as a 
      # correction is added to deal with the true centre vs the centre used to create the mesh points, where sides are of different 
      # length, the elevation won't increase immediately from the known input centre elevation as the distances from the centre are arbitrary
      # --- furthermore, the centre elevation isn;t the deepest - where one edge elevation is lower than the other, the deepest part of the 
      # parabola will actually be closer to the lower edge
      
      row_est=400
      norm_elev_matrix <- matrix(NA, nrow=row_est, ncol=1) # pre-define list size to hold values
      
      #norm_r_elev
      for (x in 1:nrow(cl_eg)){ 
        
        #x_val=cl_eg[x,6] # <<< distance from centreline (if centreline is not at 0, these values are just part of the overall x scale for which the paraobola was calculated - no need to correct them as the x axis need only be internally consitent)
        x_val=cl_eg$cl_dist_r[x]
              
        #side<-cl_eg[x,9]  # <<< side value
        side=cl_eg$side[x]

        if(verbose==1){cat("Side: ", side, "\n")}
        
        # if to the side is 1 (i.e. right) then multiply x_val by -1 (so distances run from negative to positive as they move across the centreline)
        if(side == 1){
          #if (true_centre==TRUE){
          #  x_val=(x_val*-1)#+correction
          #} else {
            x_val=(x_val*-1)
          #}
        } else if (side == 2){
          #if (true_centre==TRUE){
          #  x_val=(x_val)#+correction
          #} else {
            x_val=(x_val)
          #}
        }

        if(verbose==1){cat("x_val: ", x_val, "\n")}
        r_depth = (a*(x_val^2))+(b*x_val)+c
        
        r_depth=round(r_depth)
        if(verbose==1){cat("r_depth: ", r_depth, "\n")}
        
        norm_elev_matrix[x]<-rbind(r_depth)
      }
      
      norm_elev_matrix<-na.omit(norm_elev_matrix) # get rid of nans
      norm_elev_matrix<-matrix(norm_elev_matrix, nrow=length(norm_elev_matrix), byrow=TRUE) # reshape it
      cl_elev=(a*(x_val^2))+(b*x_val)+c
      
      ##################
      # Add data to data frame (a new row for each point) 
      for (i in 1:nrow(cl_eg)){
        
        # This step uses knowledge of the last index used by the last line to ensure data isn't overwritten
        if(count==0){
          pos=count+i   
        }
        if(count>=1){
          #if(verbose==1){print("Count greater than 1")}
          pos=last_indx_counter[count]+i    
          #pos=last_indx_counter+i    
        }
        
        if(verbose==1){print(pos)}

        out_df[pos,]<-cbind(cl_eg[i,], norm_elev_matrix[i], cl_elev)
     
        
      }
      
      count=count+1
      
      #Append last index to ensure lines of out_df in next loop iteration aren't overwritten
      if (count==1){
        last_indx_counter[count]<-rbind(nrow(cl_eg))
      }
      
      if (count>1){
        last_indx_counter[count]<-rbind(nrow(cl_eg)+last_indx_counter[count-1])
      }

      

      outputs=list("out_df"=out_df, "count"=count, "last_indx_counter"=last_indx_counter, "true_cl_from_zero"=mid_depth[1])
      return(outputs)
    } else {
      return(outputs=FALSE)
    }
  }
}

#' true_centre=TRUE # If true, the parabola is calculated for the true centre for a given normal, calculated as the point exactly half way between the fjord edges - this often differs 
#'                    from the mapped centreline
#'                    Using this flag, the output dataframe includes 3 additional columns defining the TRUE cl (x,y) and offset from the mapped cl (x,y) and also returns a plot showing 
#'                    the elevations nearest to the true centre predicted using the piece-wise approach (as well as the actual TRUE centre points used to define the piece-wise model)
#'
#' OUTPUTS
#'
#' elev_mesh = normal and centreline points with assigned elevations
#'           FORMAT IF (true_centre==TRUE): "x","y","cl_x","cl_y","cl_ID","cl_dist_r","cl_dist_l", "path", "side","norm_elev","cl_elev", "true_cl_x", "true_cl_y", "dist_from_cl"
#'           FORMAT IF (true_centre!=TRUE): "x","y","cl_x","cl_y","cl_ID","cl_dist_r","cl_dist_l", "path", "side","norm_elev","cl_elev"
#'
#' xyz       = points used to define the piece-wise model - if true_centre==TRUE, then these will be the EXACT centre (x,y) (not the nearest mesh node equvalents), otherwise they will be the mapped cl(x,y) positions
#'
#' dist_threshold = observtions within this distance are associated with a given normal
channel_parabola__piecewise <- function(normal_points_FILE, edge_elevations_FILE, observations_xyz=data.frame(), f_out, seed_xyz, mouth_xyz, verbose=0, plotting=FALSE, specific_mouth_fid='', true_centre=TRUE, piecewise=1, dist_threshold=500){

  #verbose=0 # 0 turns print statements ON | 1 turns print statements OFF  
  print("Running parabola code -- will incorporate any known observations...")
  #print("this is the right code :)")

  ####################
  # read in data 
  df_norm=read_me_df(normal_points_FILE)   # x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side
  df_edge=read_me_df(edge_elevations_FILE) # X..x,y,cl_x,cl_y,cl_id,cl_dist_r,cl_dist_l,path,side,bed.elev.m.
  
  print("Checking observation_xyz input...")
  if(is.data.frame(observations_xyz)!=TRUE){stop("observations_xyz should be a n empty dataframe or of the form x,y,z -- if you have no data, do not set this variable")}
  df_obs=observations_xyz   # if observations are provided nrow will be >1, otherwise it will be 0
  print("Observation_FILE input checked.")

  # add seed and mouth to obs .... location...
  seed_elevation=seed_xyz[3]
  mouth_elevation=mouth_xyz[3]

  if (plotting==TRUE){
    
    df_norm_temp=df_norm
    df_obs_temp=df_obs
    seed_xyz_temp=seed_xyz
    mouth_xyz_temp=mouth_xyz
    
    df_norm_temp$type="01:mesh"
    seed_xyz_temp$type="02:seed"
    mouth_xyz_temp$type="03:mouth"
    df_obs_temp$type="04:obs"
    
    ggplot(df_norm_temp, aes(x,y,colour=type))+
      geom_point() + 
    #  geom_point(data=df_obs_temp, aes(x,y,colour=type))+ 
      geom_point(data=seed_xyz_temp, aes(x,y,colour=type)) + 
      geom_point(data=mouth_xyz_temp, aes(x,y,colour=type)) +
      ggtitle("Mesh vertices and observations")
  }

  ####################
  # Clips the normal point dataframe up to a declared mouth fid 
  # If a specific mouth fid is declared, then ignore all points associated with fid values > specific_mouth_fid
  # e.g specific_mouth_fid=589 # all oints > fid 589 will be ignored
  if (specific_mouth_fid!=''){
    cat("Subsetting to cl_id: ", specific_mouth_fid, "\n")
    df_norm<-df_norm[df_norm$cl_id<=specific_mouth_fid,] 
    df_edge<-df_edge[df_edge$cl_id<=specific_mouth_fid,] 
  }


  ####################
  # get number of unique lines (based on fid)
  print("Querying number of lines...")
  query_df<-subset(df_norm, !duplicated(cl_id))
  unique_cl_norms=unique(df_norm$cl_id)
  mn_cl_id<-min(query_df$cl_id)
  mx_cl_id<-max(query_df$cl_id)

  
  # strip edge file to same cl ids as in the norm file (in case you have manually altered the norm file)
  print("Getting min and max edge IDs...")
  df_edge=df_edge[is.element(df_edge$cl_id, unique_cl_norms),]
  
  query_edge_df<-subset(df_edge, !duplicated(cl_id))
  mn_edge_cl_id<-min(query_edge_df$cl_id)
  mx_edge_cl_id<-max(query_edge_df$cl_id)


  ####################
  # calc fid interval (we want it to be 1)
  interval=(length(query_df$cl_id)-1)/(mx_cl_id-mn_cl_id) # check FID interval - check if 1 by comparing length (minus 1!) compared to diff between min and max FID values

  if(interval!=1){
    print("/n/nWARNING: centreline points do not have an FID interval of 1/n/n")
    numPoints=nrow(query_df) # the df row count is equal to the unique cl_id values, regardless of whether the cl_id values themselves are sequential
    #stop("Need to manually define loop step as centreline points do not have an FID interval of 1")
  } else{
    numPoints=(mx_cl_id-mn_cl_id)+1 # if cl_id interval is equal to 1, then use the cl_id range to define the number of centreline elevation points
  }

  # CHECK THAT EDGE ELEVATION FILE HAS THE SAME NUMBER OF CENTRELINE FIDS (OTHERWISE IT ISN'T THE RIGHT FILE) 
  #  expect_matching_fid(mx_cl_id,mn_cl_id,mx_edge_cl_id,mn_edge_cl_id)
  cl_ids=length(seq(from=mn_cl_id, to=mx_cl_id, by=1))

  #######################
  # Calculate centreline elevations based on parabola informed by edges and observations - the centreline of the parabola at each transect ID is then used to inform the linear model
  # cl_elevs_df contains a node for each centreline fid (true or mapped depending on setting of true_centre) -- if there is no observation near a normal, then the z value for the point is NA
  # This is then considered further down the line to ensure predictive models consider the correct cumulative spacing between points, handling z=NA instances accordingly.
  #print("Entering the main loop...")
  if (nrow(df_obs)>1){ # if no obs, df_obs will have 2 rows -> the seed and the mouth

    
    print("**************************")
    print("Observations available.")
    print("**************************")
    print("About to call parabola_centreline_elev_from_OBS()")
    print(paste0("Min observation: ", min(df_obs$z)))
    print(paste0("Max observation: ", max(df_obs$z)))

    cl_elevs_df=parabola_centreline_elev_from_OBS(df_norm, df_edge, df_obs, verbose=verbose, plotting=FALSE, dist_threshold=dist_threshold, true_centre=true_centre)

    print("Finished with parabola_centreline_elev_from_OBS()")

    #check how many centreline nodes could be populated using obs
    num_pred_cl_elev=nrow(cl_elevs_df)
    print(paste0("Number of centreline elevations assigned from observations: ", num_pred_cl_elev))
    
    if (plotting==TRUE){  
      cl_elevs_df_TEMP=cl_elevs_df
      df_norm_TEMP=df_norm
      df_edge_TEMP=df_edge
      cl_elevs_df_TEMP$type="True cl xy (between edges)"
      df_norm_TEMP$type="Mapped cl xy"
      df_edge_TEMP$type="Edge"
      ggplot(df_norm_TEMP, aes(cl_x,cl_y, colour=type)) +
        geom_point() +
        geom_point(data=cl_elevs_df_TEMP, aes(x,y, colour=type))+
        geom_point(data=df_edge_TEMP, aes(X..x,y, colour=type))+
        geom_text(data=cl_elevs_df_TEMP, aes(x,y,label=centre_fids),size=2,check_overlap=TRUE,vjust=0.7,hjust=-0.3,nudge_x=0, colour="black") +
        ggtitle("True centreline points vs centreline mapped nodes")
    }

    # Combine, seed, cl and head xyz into one dataframe
    cl_elevs_df_SD_FILTERED=cl_elevs_df[abs(cl_elevs_df$centre_elevs) > (sd(cl_elevs_df$centre_elevs)),]
    print(paste0("Number of centreline elevations assigned from observations AFTER filtering: ", nrow(cl_elevs_df_SD_FILTERED)))

    xyz=combine_seed_mouth_cl_elevs(cl_elevs_df_SD_FILTERED, seed_xyz, mouth_xyz, verbose=verbose, cumdist=FALSE)
    
  } else {
    
    if(verbose==1){
      print("**************************")
      print("No observations available.")
      print("**************************")
    }

    xyz=combine_SPARSE_seed_centre_mouth_xyz(seed_xyz, df_norm, mouth_xyz, df_edge, plotting=plotting, true_centre=true_centre, verbose=verbose, cumdist=FALSE)
    cl_elevs_df=xyz[2:(nrow(xyz)-1),] # all centre points excluding the seed and mouth points
                                    # the variable cl_elevs_df has a different header if obs are available but only the x and y columns are required by the other steps
   if (nrow(xyz)!=(cl_ids+2)){
     stop("xyz should be 2 rows longer than the total number of unqiue centreline nodes as\nthe seed and mouth points have now been included...\nthis is not the case...\ncheck channel_parabola__piecewise()")
   }
 }


  #threshold=20 # should be passed into the overall function at the top
  #if (num_pred_cl_elev>=(length(cl_ids)-((length(cl_ids)/100)*threshold)) | piecewise==1){
  if (piecewise==1){
    #print("Number of centreline elevations estimated using observations exceeds threshold - using piece-wise approach...")
    if(verbose==1){
      print("***************************************************************")
      print("Piece-wise linear method used to calculate centreline elvations")
      print("***************************************************************")
      }
    model_required='piecewise'
    
  } else {
    if(verbose==1){
      print("***************************************************************")
      print("Number of centreline elevations estimated using observations under threshold - using linear model...")
      print("***************************************************************")
      }
    model_required='linear'
    
    # calc coefficients for linear model
    lm_out=linear_model_elevation_coeffs(xyz, plotting=TRUE) # should take in original seed and mouth xyz?
    lm_coeffs=lm_out[[1]]
    if (plotting==TRUE){lm_out[2]}
  }

  ###################
  # Populate mesh with elevations
    
  # Prefdefine output
  # See here: http://stackoverflow.com/questions/3642535/creating-an-r-dataframe-row-by-row
  N <- 1e4  # some magic number (possibly an overestimate) for pre-setting rows of dataframe
  
  out_df=data.frame(x=rep(NA, N), y=rep(NA, N), cl_x=rep(NA, N), cl_y=rep(NA, N), cl_id=rep(NA, N), cl_dist_r=rep(NA, N), cl_dist_l=rep(NA, N), path=rep(NA, N), side=rep(NA, N), norm_elev=rep(NA, N), cl_elev=rep(NA, N), stringsAsFactors=FALSE) 
  
  last_indx_counter=matrix(NA, nrow=N, ncol=1) # an empty list used to hold the last index filled by each dataset in set_elevations()
  
  cumulative_dist_between_centrenodes=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  centre_node_x=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  centre_node_y=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  count=0

  elevs_to_check=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  true_centre_x_CHECK=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  true_centre_y_CHECK=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  dists_CHECK=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  fid_CHECK=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)
  dist_from_mapped_centre=matrix(NA, nrow=(mx_cl_id-mn_cl_id), ncol=1)

  problem_normals=data.frame(cl_id=double())


  # Assign distances to centrepoints based on the mapped centreline 
  #   - this is important as the centrepoints (from obs) may not cover the total 
  #     length of the chaannel and if you only use the cumulative dist of those points, 
  #     you will underestimate length and end up having all sorts of problems, especially 
  #     if estimating along channel centreline elevations using sparse nodes close to observations

  print("channel_parabola__piecewise(): Assigning cumulative dist from mapped centreline to nodes close to observations")
  df_norm_cl=subset(df_norm, !duplicated(df_norm$cl_id) )
  df_norm_cl=df_norm_cl[c('cl_id', 'cl_x','cl_y')]
  mapped_channel_cumulative_dist=cumulative_dist(df_norm_cl, use_cl_xy=TRUE)

  print("channel_parabola__piecewise(): cumulative_dist() complete")
  print("channel_parabola__piecewise(): Calculating nearest observation to node...")

  cum_dist_collection=c()
  for (node_with_ob in 1:nrow(xyz)){
      
      #print(node_with_ob)
      #print(xyz[node_with_ob,])
      node_with_ob_inst=xyz[node_with_ob,]

      nn2=get.knnx(as.matrix(node_with_ob_inst[c('x','y')]), as.matrix(mapped_channel_cumulative_dist[c('cl_x','cl_y')]), k=1, algorithm=c("kd_tree")) 
      min_dist_of_norm=min(nn2$nn.dist)       # value of shortest distance
      min_indx_dist_of_norm=which.min(nn2$nn.dist) # index of shortest distance
      nearest_mapped_centrenode=mapped_channel_cumulative_dist[min_indx_dist_of_norm,]

      # plot to check...
      
      # Assign to dataframe... (at correct row index!)
      #xyz$cumulative_dist= ## use nearest neighbout - only take in x and y!!
      cum_dist_collection[node_with_ob]=nearest_mapped_centrenode$cumulative_dist
  }

  xyz$cumulative_dist=cum_dist_collection  

  #plot(df_norm_cl$cl_x, df_norm_cl$cl_y, col="red")
  #points(xyz$x, xyz$y, col="blue")
  print("channel_parabola__piecewise(): Calculating nearest observation to node complete.")

  print("channel_parabola__piecewise(): Assigning elevations...")
  for (fid in mn_cl_id:mx_cl_id){

    #print(paste0("channel_parabola__piecewise(): Working on: ", fid))

    fid_CHECK[count+1]=fid

    df_norm_sub=subset(df_norm, cl_id==fid)

    if(nrow(unique(df_norm_sub))<2){
        cat("XXXXXXXXXXXX\n")
        cat("First catch for bad cl_ids...\n")
        cat("The normal along this cl_id will not have any elevations and will be omitted from the output - possible that the normal consists of < 3 points\n")
        cat("XXXXXXXXXXXX\n")

          # STILL ACCUMUALTE DISTANCE - THE POINT STILL EXISTS IN SPACE!
          POINT1_XY=data.frame('x'=centre_node_x[count],'y'=centre_node_y[count],'z'=0) # the last point  <<<< need to change this in case the last point was NA
          POINT2_XY=data.frame('x'=centre_node_x[count+1],'y'=centre_node_y[count+1],'z'=0) # the current point
          dist_along_centreline=cumulative_dist(rbind(POINT1_XY, POINT2_XY), verbose=verbose)[2,4]       
          cumulative_dist_between_centrenodes[count+1]= cumulative_dist_between_centrenodes[count]+dist_along_centreline # the last distance + the new distance i.e. it accumulates
          dists_CHECK[count+1]=dist_along_centreline

        count=count+1 

    } else {        

      if (fid<=max(df_edge$cl_id)){
        if(verbose==1){cat("Using normal cl_id...\n")}
        df_edge_sub=subset(df_edge, cl_id==fid)   
      } else {
        # Assumed that the mismatch in cl_id between norm and edge is because the normal extends out of the coast 
        # i.e. the normal ndoes furthest from the coast share the same edge point (for which the cl_id does not change)
        if(verbose==1){cat("Normal extends past edges (norm cl_id max > edge cl_id) so using edge max...\n")}
        df_edge_sub=subset(df_edge, cl_id==max(df_edge$cl_id))
      } 
      
      edge_r=subset(df_edge_sub, side==1)
      edge_l=subset(df_edge_sub, side==2)

        # Calculate elevation at centre of parabola
        if (true_centre==TRUE){

          # centre has (x,y) equal to the TRUE centreline node

          true_centre_xy=find_point_btwixt(edge_r, edge_l, plotting=F)
           
          # TEST :: check cl_xy match those in model 
          true_centre_x_CHECK[count+1]=true_centre_xy[[1]]
          true_centre_y_CHECK[count+1]=true_centre_xy[[2]]

        # DEFUNCT - EDGES DON;T ALWAYS MATCH AS NORMALS SOMETIMES EXTEND FURTHER
        #  check_centre_xy_match_those_for_model(cl_elevs_df, true_centre_x_CHECK, true_centre_y_CHECK, count=count+1, verbose=verbose)  # func in parabola_funcs_where_obs_availble.r
   
          centre_node_x[count+1]=true_centre_xy[[1]]
          centre_node_y[count+1]=true_centre_xy[[2]]

        } else if (true_centre==FALSE){

          # centre has (x,y) equal to the MAPPED centreline node

          centre_node_x[count+1]=df_norm_sub$cl_x[1]
          centre_node_y[count+1]=df_norm_sub$cl_y[1]

        # DEFUNCT - EDGES DON;T ALWAYS MATCH AS NORMALS SOMETIMES EXTEND FURTHER
        #  check_centre_xy_match_those_for_model(cl_elevs_df, centre_node_x, centre_node_y, count=count+1, verbose=verbose)  # func in parabola_funcs_where_obs_availble.r
        
        }


        if (count == 0){ # i.e this is the first centreline fid being considered...
                  
          # Distance of first cl (x,y) from the seed
     
          #POINT1_XYZ=data.frame('x'=centre_node_x[count+1],'y'=centre_node_y[count+1], 'z'=0)
          #seed_xyz_temp<-seed_xyz
          #colnames(seed_xyz_temp)<-c('x','y','z')
          #temp=rbind(seed_xyz_temp, POINT1_XYZ)
          ##cumulative_dist_between_centrenodes[count+1]=cumulative_dist(temp, verbose=verbose)[2,4] # output is both points (2 rows)  - you want the dist of the second row from the 1st (held in second row, fourth column)
          #dists_CHECK[count+1]=cumulative_dist(temp, verbose=verbose)[2,4]

          cumulative_dist_between_centrenodes[count+1]=0 # start accumulating distance from here
          dists_CHECK[count+1]=0

        } else {
          
          POINT1_XY=data.frame('x'=centre_node_x[count],'y'=centre_node_y[count],'z'=0) # the last point  <<<< need to change this in case the last point was NA
          POINT2_XY=data.frame('x'=centre_node_x[count+1],'y'=centre_node_y[count+1],'z'=0) # the current point
          dist_along_centreline=cumulative_dist(rbind(POINT1_XY, POINT2_XY), verbose=verbose)[2,4]
          
          cumulative_dist_between_centrenodes[count+1]= cumulative_dist_between_centrenodes[count]+dist_along_centreline # the last distance + the new distance i.e. it accumulates
          dists_CHECK[count+1]=dist_along_centreline

        }

        # predict the centreline elevation based on distance along the centreline - if linear model is a poor fit, then this will not necessarily match the centreline values
        # calculated above using parabola_centreline_elev_from_OBS()
        
        if (model_required=='linear'){
          true_cl_elev=linear_model_elevations(cumulative_dist_between_centrenodes[count+1], as.double(lm_coeffs$coefficients[1]), as.double(lm_coeffs$coefficients[2]))
        } else if (model_required=='piecewise'){
          if(verbose==1){print("Using piece-wise")}
          #true_cl_elev=approx(xyz$cumulative_dist, xyz$z, xout=cumulative_dist_between_centrenodes[count+1], method="linear")$y # works as expected even if xyz$z contains NA values
          true_cl_elev=piecewise_model_elevation(xyz, xout=cumulative_dist_between_centrenodes[count+1])
        }
        
        elevs_to_check[count+1]=true_cl_elev # cross check with y=mx+c and cumulative_dist_between_centrenodes
      
        output=set_elevations(df=df_norm_sub, 
                           df_edge=df_edge_sub, 
                           fid=fid, 
                           out_df=out_df, 
                           count=count, 
                           centreline_depth=true_cl_elev, 
                           true_centre=true_centre,
                           last_indx_counter=last_indx_counter,
                           verbose=0)

        if (class(output)=="list"){
          out_df=output$out_df
          count=output$count
          last_indx_counter=output$last_indx_counter   
          #print("Problem NOT reached")
          dist_from_mapped_centre[count]=output$true_cl_from_zero   

        } else {

          print("XXXXXXXXXXXX")
          print("Second catch for bad cl_ids...")
          print("The normal along this cl_id will not have any elevations and will be omitted from the output - possible that the normal consists of < 3 points")
          print("XXXXXXXXXXXX")
          
          # No data should be added for this FID - equally, it is still a node so the accumualted dist need not change
          problem_normals[nrow(problem_normals)+1,]=fid

            # STILL ACCUMUALTE DISTANCE - THE POINT STILL EXISTS IN SPACE!
            POINT1_XY=data.frame('x'=centre_node_x[count],'y'=centre_node_y[count],'z'=0) # the last point  <<<< need to change this in case the last point was NA
            POINT2_XY=data.frame('x'=centre_node_x[count+1],'y'=centre_node_y[count+1],'z'=0) # the current point
            dist_along_centreline=cumulative_dist(rbind(POINT1_XY, POINT2_XY), verbose=verbose)[2,4]       
            cumulative_dist_between_centrenodes[count+1]= cumulative_dist_between_centrenodes[count]+dist_along_centreline # the last distance + the new distance i.e. it accumulates
            dists_CHECK[count+1]=dist_along_centreline
          
          count=count+1 # required to ensure cumualtive distance keeps increasing - the 'problem' cl_id still exists, just the elevations are going to be ignored


        }
    }
  }
  print("channel_parabola__piecewise(): Assigning elevations complete.")

  ##################
  ## Output

    # check true elevs are held within the output - these will not match perfectly as:
      # although the true centre is important for defining the [parabola but the specific point between the bank edges is not likely cleanly 
      # divisbile by 2 - thus, as the bank distances from 0 have been used to calc the parabola (setting its x xaxis range), we now 
      # have to calculate a value at 0 to have a complete chain of points - whether it is the centre or not is now irrelevant
      # -- see x_val setting for defing cl_elev in set_elevations() (to have the following plot straight, this x_val would be assigned the 
      # precise midpoint between the two edges which does not necessaily match the spacing of the points along the normal)
     

      # only makes sense if the "0" cl was used
      #output_cl_elevs<-out_df[!duplicated(out_df$cl_id),]$cl_elev # cl_elevs within the output dataframe - these are calculated using the parabola function for the mid point of each normal
      #plot(elevs_to_check-output_cl_elevs, main="May not match perfectly as centreline elevs do\n not account for normal node spacing", x_lab="Input centreline elev", y_lab="Output centreline elev")

  out_df<-out_df[!duplicated(out_df),]
  out_df<-na.omit(out_df) 

  write.csv(out_df, file=f_out, row.names=FALSE)

  # Write out path fids that had issues
  if (nrow(problem_normals)>0){
    f_out_ISSUES=paste0(file_path_sans_ext(f_out), "___CL_ID_ISSUES.csv")
    write.csv(out_df, file=f_out_ISSUES, row.names=FALSE)    
  }

  # Add cols defining true_cl_x and true_cl_y + dist of true from mapped cl to out_df
  # These values are used to check the elvations predicted in set_elevations - see notes 
  # at top of check_obs_using_piecewise()
  if (true_centre==TRUE){

    fid_CHECK<-na.omit(fid_CHECK) # list of fids
    true_centre_x_CHECK<-na.omit(true_centre_x_CHECK) # list of true cl x for each fid
    true_centre_y_CHECK<-na.omit(true_centre_y_CHECK) # list of true cl y for each fid
    id_true_cl_xy_df=data.frame("fid"=matrix(fid_CHECK),"true_cl_x"=matrix(true_centre_x_CHECK),"true_cl_y"=matrix(true_centre_y_CHECK), "dist_from_mapped_centre"=matrix(dist_from_mapped_centre))
    
    count=0
    for (fid_inst in mn_cl_id:mx_cl_id){
      count=count+1
      out_df$true_cl_x[out_df$cl_id==fid_inst]<-id_true_cl_xy_df$true_cl_x[id_true_cl_xy_df$fid==fid_inst]
      out_df$true_cl_y[out_df$cl_id==fid_inst]<-id_true_cl_xy_df$true_cl_y[id_true_cl_xy_df$fid==fid_inst]
      out_df$dist_from_cl[out_df$cl_id==fid_inst]<-id_true_cl_xy_df$dist_from_mapped_centre[id_true_cl_xy_df$fid==fid_inst]
    }

    colnames(out_df)<-c("x","y","cl_x","cl_y","cl_id","cl_dist_r","cl_dist_l", "path", "side","norm_elev","cl_elev", "true_cl_x", "true_cl_y", "dist_from_cl")  
    
    #sanity_check_plot=check_obs_using_piecewise(xyz, out_df, seed_xyz, mouth_xyz) # see parabola_funcs_where_obs_availble.r

    return(list(out_df, xyz))    

  } else {
   
    return(list(out_df, xyz))    
  
  } 
   
}

