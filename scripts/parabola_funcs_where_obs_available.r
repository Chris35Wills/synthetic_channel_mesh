library(FNN)
library(sp)
library(ggplot2)
library(raster)
library(foreach)

read_me_df<-function(data_f){
  data=read.csv(data_f)
  df<-as.data.frame(data)
  return(df)
}

# edge_r/edge_l - dataframe containg x and y of each edge point (col names passed using edge_x and edge_y)
# centre - dataframe containg x and y of centre point (cols must be 'x' and 'y' -- can;t get colour setting to work using aes_string...)
plot_edge_and_mid_points<-function(edge_r, edge_l, centre, edge_x='x', edge_y='y', centre_x='x', centre_y='y') {
  ggplot(edge_r, aes_string(edge_x, edge_y)) +
        geom_point() +
        geom_point(data=edge_l, aes_string(edge_x, edge_y))+
        geom_point(data=centre, aes(x, y, col='red'))
}

# Get coordinate of point that lies between 2 known points
# pnt_1_xy - of type list where x and y coords are held in pnt_1_xy[1] and pnt_1_xy[2] respectively
# pnt_2_xy - as for pnt_1_xy
find_point_btwixt<-function(pnt_1_xy, pnt_2_xy, plotting=TRUE){

  true_centre_xy_diff=c(((pnt_1_xy[1]-pnt_2_xy[1])/2),((pnt_1_xy[2]-pnt_2_xy[2])/2))
  true_centre_xy=c((pnt_1_xy[1]-true_centre_xy_diff[1]),(pnt_1_xy[2]-true_centre_xy_diff[2]))
  #true_centre_xy=c((pnt_1_xy[1]-true_centre_xy_diff[1]),(pnt_1_xy[2]-true_centre_xy_diff[2]))

  if (plotting==TRUE){  
    plot(-8:8,-8:8,type="n")
    points(pnt_2_xy[1],pnt_2_xy[2], col="blue", pch=19)
    points(pnt_1_xy[1],pnt_1_xy[2], col="red", pch=19)
    points(true_centre_xy[1],true_centre_xy[2], col="green", pch=17)
    text(pnt_1_xy[1],pnt_1_xy[2], labels="point 1", cex= 0.7, pos=3)
    text(pnt_2_xy[1],pnt_2_xy[2], labels="point 2", cex= 0.7, pos=3)
    text(true_centre_xy[1],true_centre_xy[2], labels="btwixt", cex= 0.7, pos=3)
  }

  return(true_centre_xy)
}


#' Reads in the normal and edge data
#' Reads in some observations
#' Loops through each centreline and finds the nearest observation to the nodes that share the same centreline ID
#' Where an obserrvation is within a threshold_distance of any points along the normal, this is then matched to the nearest point along the normal (where the normal point adopts the elevation of the observation - we call this the pseudo-observation)
#' Using the known edge elevations and the pseudo-elevation, the parabola variables are calculated
#' The elevation at the centre of the normal is then claculated (i.e. half way between the 2 edge nodes as opposed to the "0" centreline distance within the normal point dataframe, as this is from the digitisation process) 
#' The centre elevation and the cl_fid whith which it is associated is then stored
#' 
#' These centereline elevations can then be used to calculate a model for the down fjord elevation profile (also considering a mouth and seed elevation)
#' The mesh is then populated as using channel_parabola_slope() except rather than using a linear trend between mouth and seed elevations, a real trend can now be provided
parabola_centreline_elev_from_OBS <- function(df_norm, df_edge, df_obs, verbose=0, plotting=FALSE, dist_threshold=2000, true_centre=TRUE, min_expected=-6000, max_expected=0){

  #verbose=0 # 0 turns print statements ON | 1 turns print statements OFF  
  if(verbose==1){cat("Running parabola code where observations are known...\n")}

  # Calculate centreline ID range
  cl_id_min=min(df_norm$cl_id)
  cl_id_max=max(df_norm$cl_id)

  centre_elevs=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
  centre_fids=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
  centre_x=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
  centre_y=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)

  # subset df_norm to df_edge       
  count=0
  for (cl_id in cl_id_min:cl_id_max){
    
    if(verbose==1){cat("parabola_centreline_elev_from_OBS() -- working on cl_id: ", cl_id, "\n")}
    
    # subset points
    norm_inst=df_norm[df_norm$cl_id==cl_id,]
    
    if (cl_id<=max(df_edge$cl_id)){
      if(verbose==1){cat("Using normal cl_id...\n")}
      edge_inst=df_edge[df_edge$cl_id==cl_id,]
    } else {
      # Assumed that the mismatch in cl_id between norm and edge is because the normal extends out of the coast 
      # i.e. the normal ndoes furthest from the coast share the same edge point (for which the cl_id does not change)
      if(verbose==1){cat("Normal extends past edges (norm cl_id max > edge cl_id) so using edge max...\n")}
      edge_inst=df_edge[df_edge$cl_id==max(df_edge$cl_id),]
    }


  if (nrow(edge_inst)<2){
    cat("\nxxxxXXXXxxxx\n")
    cat("An edge point has been lost -- possibly when masking\n")
    cat("This only affects cl_id ", cl_id, "\n")
    cat("This normal will be skipped\n")
    cat("xxxxXXXXxxxx\n\n")
  } else {

    if (min(norm_inst$cl_id)!= cl_id | max(norm_inst$cl_id) != cl_id) {print("We have a problem - not selecting correct normal based on cl_id")}

    if (round(unique(norm_inst$cl_x)) != round(unique(edge_inst$cl_x)) | round(unique(norm_inst$cl_y)) != round(unique(edge_inst$cl_y))) {
      cat("norm (cl_x,cl_y): (", norm_inst$cl_x[1],",",norm_inst$cl_y[1], ")\n", sep="")
      cat("edge (cl_x,cl_y): (", edge_inst$cl_x[1],",",edge_inst$cl_y[1], ")\n", sep="")
      stop("We have a problem - edge and norm cl_xy don't match - see parabola_centreline_elev_from_OBS()")
    }

    # get nearest observation to mapped centre node 
    #keeps=c('x', 'y', 'cl_id')
    #print("Only considering nearest observation to centre of the current normal...")
    #print("Getting nearest ob to centrenode...")

    keeps=c('cl_x', 'cl_y', 'cl_id') # <<<<<<<<<<<<< JUST CONSIDER CENTRE
    df_xy=norm_inst[keeps]    # << trim norm_inst to just three columns for the purposes of running the nn analysis >> df_xy has same rows a norm_inst
    df_xy<-df_xy[!duplicated(df_xy),]

    
      # get.knnx must take ONLY the xy columns
      df_xy2=df_xy
      colnames(df_xy2)<-c('x','y','z')
      df_xy2=df_xy2[c('x','y')]
      df_obs2=df_obs[c('x','y')]
    
    nn=get.knnx(as.matrix(df_xy2), as.matrix(df_obs2), k=1, algorithm=c("kd_tree")) # output will be a list the same length as df_obs
    min_dist=min(nn$nn.dist)       # value of shortest distance
    min_indx=which.min(nn$nn.dist) # index of shortest distance
    
    if (min_dist<=dist_threshold){

      ###################
      # Get nearest observation to centre
      index_nearest_ob_to_centre=nn$nn.index[min_indx] # index of normal point closest to nearest observation
      nearest_ob_to_centre=df_obs[min_indx,]                                   # xyz of nearest observation
     
      if (is.na(nearest_ob_to_centre$z)==FALSE & (nearest_ob_to_centre$z < max_expected) & (nearest_ob_to_centre$z > min_expected)){ # avoid any spurious values

        ###################
        ## Get nearest normal to observation
        keeps=c('x','y','side')
        df_norm_xy=norm_inst[keeps]
        
        #nn2=get.knnx(nearest_ob_to_centre, df_norm_xy, k=1, algorithm=c("kd_tree")) 
        nn2=get.knnx(as.matrix(nearest_ob_to_centre[c('x','y')]), as.matrix(df_norm_xy[c('x','y')]), k=1, algorithm=c("kd_tree")) 
        min_dist_of_norm=min(nn2$nn.dist)       # value of shortest distance
        min_indx_dist_of_norm=which.min(nn2$nn.dist) # index of shortest distance
        nearest_normal_point_to_ob=norm_inst[min_indx_dist_of_norm,]

        ###################
        ## Create the pseudo ob (normal node whoch has adopted the observation elevation)
        pseudo_ob=c(nearest_normal_point_to_ob$cl_dist_r,     # dist from centre line of nearest point along the normal to the observation)
                    nearest_ob_to_centre$z,                   # elev of the observation nearest the point on the normal
                    nearest_normal_point_to_ob$side)          # side of nearest point along the normal to the observation)

        # if the side of psuedo ob is 1 (i.e. right), then distance from centreline is negative (i.e. *-1) which matches the convention used throughout - postivie distance from centreline on left, negative on the right
        if (pseudo_ob[3]==1){
          pseudo_ob[1]=pseudo_ob[1]*-1
        } 

        ###################
        # Define right and left edges of parabola
        right_max_cl_dist=max(edge_inst[edge_inst$side==1,]$cl_dist_r)
        right_edge_elev=subset(edge_inst$bed.elev.m., edge_inst$side==1)
        left_max_cl_dist=max(edge_inst[edge_inst$side==2,]$cl_dist_l)
        left_edge_elev=subset(edge_inst$bed.elev.m., edge_inst$side==2)
             
        edge_r=c(right_max_cl_dist*-1,right_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #1
        edge_l=c(left_max_cl_dist,left_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #2

        if (as.integer(pseudo_ob[1]) != as.integer(edge_r[1]) & as.integer(pseudo_ob[1]) != as.integer(edge_l[1])){  # this happens if an observation is on the edge of a channel - we ignore these points (and use the GIMP DEM instead)

          if(verbose==1){cat("Right edge elevation: ", right_edge_elev, "m\n", "Left edge elevation: ", left_edge_elev, "m\n", sep="")}
            
          #x will be distance and y will be elevation...
          abc=calc_parabola_vertex(edge_r[1], edge_r[2], pseudo_ob[1], pseudo_ob[2], edge_l[1], edge_l[2])
          
          ###################
          # Extract centre point (true or mapped) elevation from parabola 
          if (true_centre==FALSE){
            
            ############
            # calculate just the "0 m from cl" elevation (this is the node on the digitized centreline)
            digitised_centre_dist=0
            digitized_centre_depth=round(abc[1]*(digitised_centre_dist^2))+(abc[2]*digitised_centre_dist)+abc[3]
            
            # append cl_id and elev to list
            count=count+1  
            centre_elevs[count]=digitized_centre_depth
            centre_x[count]=norm_inst$cl_x[1]
            centre_y[count]=norm_inst$cl_y[1]
            centre_fids[count]=cl_id

            # print(paste0("cl_id nearest to observation: ", cl_id))
            # print(paste0("Nearest observation: ", digitized_centre_depth))

          } else if (true_centre==TRUE){

            ############
            # calculate the true centre elevation (this is the centre point relative to the two known edge points and will likely differ from the digitized centreline)
            true_centre_dist=(edge_l[1]+(edge_r[1]))/2 # in units of distance from the digitized centre
            true_centre_depth=round(abc[1]*(true_centre_dist^2))+(abc[2]*true_centre_dist)+abc[3]
            
            # calculate the xy coordinate of the centre point
            edge_l_xy=c(edge_inst[edge_inst$side==2,]$X..x, edge_inst[edge_inst$side==2,]$y)
            edge_r_xy=c(edge_inst[edge_inst$side==1,]$X..x, edge_inst[edge_inst$side==1,]$y)
            true_centre_xy=find_point_btwixt(edge_l_xy, edge_r_xy, plotting=plotting)

            # append cl_id and elev to list
            count=count+1  
            centre_elevs[count]=true_centre_depth
            centre_x[count]=true_centre_xy[1]
            centre_y[count]=true_centre_xy[2]
            centre_fids[count]=cl_id

            # print(paste0("Nearest observation: ", true_centre_depth))

          }
        }

      }
    } else {
     
      #cat("Nearest observation to normal is > than the distance threshold set: ", dist_threshold, "\n", sep="")
      if(verbose==1){cat("Nearest observation to normal is > than the distance threshold set: ", dist_threshold, "\n", sep="")}
      if(verbose==1){cat("Moving to next normal.\n", sep="")}

      # still need to write the data out - z is equivalent of NA (set here to -9999)
      if (true_centre==FALSE){
        
        ############
        # calculate just the "0 m from cl" elevation (this is the node on the digitized centreline)
        digitised_centre_dist=0
                
        # append cl_id and elev to list
        count=count+1  
        centre_elevs[count]=-9999
        centre_x[count]=norm_inst$cl_x[1]
        centre_y[count]=norm_inst$cl_y[1]
        centre_fids[count]=cl_id

      } else if (true_centre==TRUE){

        ############
        # calculate the true centre elevation (this is the centre point relative to the two known edge points and will likely differ from the digitized centreline)
        true_centre_dist=(edge_l[1]+(edge_r[1]))/2 # in units of distance from the digitized centre
        
        # calculate the xy coordinate of the centre point
        edge_l_xy=c(edge_inst[edge_inst$side==2,]$X..x, edge_inst[edge_inst$side==2,]$y)
        edge_r_xy=c(edge_inst[edge_inst$side==1,]$X..x, edge_inst[edge_inst$side==1,]$y)
        true_centre_xy=find_point_btwixt(edge_l_xy, edge_r_xy, plotting=plotting)

        # append cl_id and elev to list
        count=count+1  
        centre_elevs[count]=-9999
        centre_x[count]=true_centre_xy[1]
        centre_y[count]=true_centre_xy[2]
        centre_fids[count]=cl_id

      }

    }

  }

  }
#  centre_elevs=na.omit(centre_elevs)#[,1]
#  centre_fids=na.omit(centre_fids)#[,1]
#  centre_x=na.omit(centre_x)#[,1]
#  centre_y=na.omit(centre_y)#[,1]

  #centre_elevs <- within(centre_elevs, centre_elevs[centre_elevs==4] <- NA)
  if (verbose==1){cat("Writing elevs to dataframe\n")}
  elevs_fid_df=data.frame("x" = c(centre_x),
                          "y" = c(centre_y),
                          "centre_elevs"=c(centre_elevs), 
                          "centre_fids"=c(centre_fids)
                          )

  elevs_fid_df$centre_elevs[elevs_fid_df$centre_elevs==-9999] <- NA 
  # set any -9999 values to NA (these are the centre nodes with no enarest observation)
  elevs_fid_df=na.omit(elevs_fid_df)
  #ggplot(elevs_fid_df, aes(centre_fids, centre_elevs)) +
  #      geom_point() +
  #      geom_smooth(method = "lm", se = FALSE) + 
  #      ggtitle("Elevations from parabolas constructed \nwhere normals are near to observations")

  return(elevs_fid_df)
  
}

# TEST :: check cl_xy match those used for definition of piecewise model 
check_centre_xy_match_those_for_model<-function(cl_elevs_df, true_centre_x_CHECK, true_centre_y_CHECK, count, verbose=1){

  cl_xy_in_loop=data.frame("cl_x"=matrix(true_centre_x_CHECK), "cl_y"=matrix(true_centre_y_CHECK))

  if ((cl_elevs_df$x[count]!=cl_xy_in_loop$cl_x[count])|
     (cl_elevs_df$y[count]!=cl_xy_in_loop$cl_y[count])){
    print("xy dont match")
  } else {

    if(verbose==1){
      cat("\nIndex ", count, "\n", sep="")
      cat("cl_elev_x: ", cl_elevs_df$x[count], "\n", sep="")
      cat("cl_loop_x: ", cl_xy_in_loop$cl_x[count], "\n", sep="")
      cat("cl_elev_y: ", cl_elevs_df$y[count], "\n", sep="") 
      cat("cl_loop_y: ", cl_xy_in_loop$cl_y[count], "\n", sep="")
    }
  }
}

#' Calculate euclidea distance between 2 points
#' x1 or x2 can take the form c(0,0)
#' e.g. euc.dist(c(0,0), c(3,3))
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#' Calculates cumulative distance of a series of ordered x,y points
#' @param df - dataframe of format: x,y,centre_elevs,n....
#' @param use_cl_xy - if TRUE, cols used are named 'cl_x' and 'cl_y' otherwise 'x' and 'y'
cumulative_dist<-function(df, use_cl_xy=FALSE,verbose=1){
  if(verbose==1){print("Running calc_cumu_dist()...")}

  #coordinates(df) <- c("x", "y")
  #df$cumulative_dist <- spDistsN1(df, df[1,]) # df[1,] relates to the z variable
  #df=as.data.frame(df)
  
  if (use_cl_xy==FALSE){
    dists=foreach(i = 1:nrow(df), .combine = c ) %do% euc.dist(c(df[i,]$x,df[i,]$y), c(df[i-1,]$x,df[i-1,]$y))
    } else {
    dists=foreach(i = 1:nrow(df), .combine = c ) %do% euc.dist(c(df[i,]$cl_x,df[i,]$cl_y), c(df[i-1,]$cl_x,df[i-1,]$cl_y))
    }

  df$cumulative_dist=cumsum(dists)

  return(df)
}

# Combines the seed, mouth and centreline elevations
# elevs_fid_df format=('x','y','centre_elevs')
# seed_xyz/mouth_xyz format=('x','y','z')
combine_seed_mouth_cl_elevs<-function(elevs_fid_df, seed_xyz, mouth_xyz, plotting=TRUE, verbose=1, cumdist=TRUE){

    if(verbose==1){print("Running combine_seed_mouth_cl_elevs()...")}
    # combine seed, centreline elevs and mouth into one dataframe
    keeps=c('x','y','centre_elevs')
    cl_elevs=elevs_fid_df[keeps]
    keeps=c('x','y','z')
    seed_xyz=seed_xyz[keeps]
    mouth_xyz=mouth_xyz[keeps]
    colnames(cl_elevs)<-c("x","y","z")
    colnames(seed_xyz)<-c("x","y","z")
    colnames(mouth_xyz)<-c("x","y","z")
    xyz=rbind(seed_xyz, cl_elevs, mouth_xyz)
    
    if (plotting==TRUE){
      ggplot(cl_elevs, aes(x, y)) +
            geom_point() +
            geom_point(data=seed_xyz, aes(x, y, colour="red", label="seed")) +
            geom_point(data=mouth_xyz, aes(x, y, colour="green", label="mouth")) +
            ggtitle("Centreline locations") +
            coord_fixed()
    }


    ## calc distance between each neighbouring point
    if(cumdist==TRUE){
      print("*****************")
      print("WARNING: calculating cumulative dist between centrepoints - if not all centrepoints are matched or points are not equally spaced, the total cumulative distanc emay be shorted than you expect (i.e. it's 'as the crow flies')")
      print("*****************")
      xyz=cumulative_dist(xyz)
    }

    return(xyz)
}



# Combines the seed, mouth and centreline elevations - only the seed and mouth have an elevation
# Calculates the cumulative distance between these points and outputs them as x,y,z,cumulative_distance
# elevs_fid_df format=('x','y','centre_elevs')
# seed_xyz/mouth_xyz format=('x','y','z')
combine_SPARSE_seed_centre_mouth_xyz<-function(seed_xyz, df_norm, mouth_xyz, df_edge, plotting=TRUE, true_centre=FALSE, verbose=1, cumdist=TRUE){

    if(verbose==1){print("Running combine_seed_mouth_cl_elevs()...")}
    # combine seed, centreline elevs and mouth into one dataframe
    keeps=c('x','y','z')
    seed_xyz=seed_xyz[keeps]
    mouth_xyz=mouth_xyz[keeps]
    colnames(seed_xyz)<-c("x","y","z")
    colnames(mouth_xyz)<-c("x","y","z")
    

    # Calculate centreline ID range
    cl_id_min=min(df_norm$cl_id)
    cl_id_max=max(df_norm$cl_id)

    centre_elevs=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
    centre_fids=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
    centre_x=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
    centre_y=matrix(NA, nrow=(cl_id_max-cl_id_min), ncol=1)
         
    count=0
    for (cl_id in cl_id_min:cl_id_max){

            norm_inst=df_norm[df_norm$cl_id==cl_id,]
            
            # still need to write the data out - z is equivalent of NA (set here to -9999)
            if (true_centre==FALSE){

              if(verbose==1){print("Extracting mapped centre (x,y)")}              
              ############
              # calculate just the "0 m from cl" elevation (this is the node on the digitized centreline)
              digitised_centre_dist=0
                      
              # append cl_id and elev to list
              count=count+1  
              centre_elevs[count]=-9999
              centre_x[count]=norm_inst$cl_x[1]
              centre_y[count]=norm_inst$cl_y[1]
              centre_fids[count]=cl_id

            } else if (true_centre==TRUE){

              if(verbose==1){print("Extracting true centre (x,y)")}              
              edge_inst=df_edge[df_edge$cl_id==cl_id,]

              ###################
              # Define right and left edges of parabola
              right_max_cl_dist=max(edge_inst[edge_inst$side==1,]$cl_dist_r)
              right_edge_elev=subset(edge_inst$bed.elev.m., edge_inst$side==1)
              left_max_cl_dist=max(edge_inst[edge_inst$side==2,]$cl_dist_l)
              left_edge_elev=subset(edge_inst$bed.elev.m., edge_inst$side==2)
                   
              edge_r=c(right_max_cl_dist*-1,right_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #1
              edge_l=c(left_max_cl_dist,left_edge_elev) # x is the distance from centre and y is assumed zero depth < known point #2
              if(verbose==1){cat("Right edge elevation: ", right_edge_elev, "m\n", "Left edge elevation: ", left_edge_elev, "m\n", sep="")}
          

              ############
              # calculate the true centre elevation (this is the centre point relative to the two known edge points and will likely differ from the digitized centreline)
              true_centre_dist=(edge_l[1]+(edge_r[1]))/2 # in units of distance from the digitized centre
              
              # calculate the xy coordinate of the centre point
              edge_l_xy=c(edge_inst[edge_inst$side==2,]$X..x, edge_inst[edge_inst$side==2,]$y)
              edge_r_xy=c(edge_inst[edge_inst$side==1,]$X..x, edge_inst[edge_inst$side==1,]$y)
              true_centre_xy=find_point_btwixt(edge_l_xy, edge_r_xy, plotting=plotting)

              # append cl_id and elev to list
              count=count+1  
              centre_elevs[count]=-9999
              centre_x[count]=true_centre_xy[1]
              centre_y[count]=true_centre_xy[2]
              centre_fids[count]=cl_id

            }
    }

    centre_elevs[centre_elevs==-9999]=NA
    centre_xyz=data.frame("x"=centre_x,"y"=centre_y,"z"=centre_elevs)

    # combine seed, centre point and mouth
    xyz=rbind(seed_xyz, centre_xyz, mouth_xyz)
    
    ## calc distance between each neighbouring point
    if(cumdist==TRUE){
      print("*****************")
      print("WARNING: calculating cumulative dist between centrepoints - if not all centrepoints are matched or points are not equally spaced, the total cumulative distanc emay be shorted than you expect (i.e. it's 'as the crow flies')")
      print("*****************")
      xyz=cumulative_dist(xyz)
    }


    return(xyz)
}

# Test combine_SPARSE_seed_centre_mouth_xyz() produces true and mapped centreline (x,y)
# Check that only seed and mouth points have elevations
# sparse - 0,1,2
TEST_combine_SPARSE_seed_centre_mouth_xyz<-function(data_path="C:/GitHub/synthetic_channels/meshing/test_data/fenty_aoi_region/", sparse=0, verbose=0){

  dat=open_test_data(data_path, as_df=FALSE, sparse=sparse)
  df_norm=read_me_df(dat[[1]])
  df_edge=read_me_df(dat[[2]])
  df_obs=read_me_df(dat[[3]])

  seed_xyz=data.frame(x=-717233,y=-1326743,z=-158)
  mouth_xyz=data.frame(x=-700582,y=-1358375,z=-920)
  
  xyz_true=combine_SPARSE_seed_centre_mouth_xyz(seed_xyz, df_norm, mouth_xyz, df_edge, plotting=TRUE, true_centre=TRUE, verbose=verbose)
  xyz_mapped=combine_SPARSE_seed_centre_mouth_xyz(seed_xyz, df_norm, mouth_xyz, df_edge, plotting=TRUE, true_centre=FALSE, verbose=verbose)

  xyz_true$type="True"
  xyz_mapped$type="Mapped"

  ggplot(xyz_true, aes(x,y,colour=type))+
    geom_point() +
    geom_point(data=xyz_mapped, aes(x,y,colour=type)) +
    ggtitle("True centre nodes will likely be more wavey")


  x11()
  ggplot(xyz_true, aes(cumulative_dist, z, colour=type))+
    geom_point(size=5)+
  geom_point(data=xyz_mapped, aes(cumulative_dist, z, colour=type)) +
  ggtitle("Only the seed and mouth should have elevations")
}

#' Takes in data frame of form: x y z cumulative_dist (see combine_seed_mouth_cl_elevs)
#' Calculates the cummulative distance between points
#' Calcs linear model: elev~cumulative dist
#' Returns the linear model coefficients + a plot of the fit relative to the data
linear_model_elevation_coeffs<-function(xyz, plotting=TRUE){

  # calc linear model - include seed, centreline and mouth elevations
  if (any(is.na(xyz))==TRUE){print("xyz$z contains some NA values - position will still be considered in model")}

  lm_plot=ggplot(xyz, aes(cumulative_dist, z)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE)  +
        ggtitle("elevation vs cumlative point distance")

  # linear model
  #coef=lm(xyz$z~xyz$cumulative_dist)
  coef=lm(xyz$z~xyz$cumulative_dist, na.action=na.exclude) # deals with missing values should any of xyz contain z=NA - see: http://stats.stackexchange.com/questions/11000/how-does-r-handle-missing-values-in-lm
  m=coef$coefficients[2]
  c=coef$coefficients[1]

  return(list(coef,lm_plot))
}


linear_model_elevations<-function(x, c, m){
    
  lin_mod<-function(x){
    y=(m*x)+c
    cat(m, "\n")
    cat(x, "\n")
    cat(c, "\n")
    return(y)
  }

  modelled_centre_elev=lin_mod(x)
  #modelled_centre_elevs=sapply(x, lin_mod) 

  return(modelled_centre_elev)

}

# Test use of approx() when passing in NA values
test_approx_dealing_with_NA<-function(){

  # Create input data - contains some NA vals for z
  xyz=data.frame("z"=c(0,-100,-20,-200,-300,-500,NA,-900), "cum_dist"=c(0,100,200,300,400,500,600,700))

  # Create model output using xyz
  model=approx(xyz$cum_dist, xyz$z, method="linear")

  # Use model to calculate elevation at inout distances
  dists=seq(from=0,to=700,by=100)
  vals=approx(xyz$cum_dist, xyz$z, xout=dists, method="linear")$y

  # Plot results
  plot(dists,vals, col="red", main="Predicted (red) using model (blue)")
  lines(model$x, model$y, col="blue")

}

piecewise_model_elevation<-function(xyz, xout, plotting=TRUE){
  
  if (any(is.na(xyz))==TRUE){print("xyz$z contains some NA values - position will still be considered in model")}
  
  # Interpolate

  #https://stat.ethz.ch/R-manual/R-devel/library/stats/html/approxfun.html
  #piece-wise linear?
  #pred=approx(xyz$cumulative_dist, xyz$z, xout=xout, method="linear")

  x <- xyz$cumulative_dist
  y <- xyz$z
  f <- approxfun(x, y, rule=2) # rule set to 2 means that predictions out of range get set to last maximum

  if(xout>max(x)){
    print("Inside piecewise_model_elevation()")
    print("Predicting out of cumulative distance range - setting depth to maximum based on cumulative dist and centreline elevations")
    print(paste0("Max cumualtive dist: ", max(x)))
    print(paste0("Provided distance: ", xout))
    }


 # par(mfrow = c(1,1))
 # plot(x, y, main = "approx(.) and approxfun(.)")
 # points(approx(x, y), col = 2, pch = "*")
 # points(approx(x, y, method = "constant"), col = 4, pch = "*")

  
 # mod=curve(f(x-10), min(xyz$cumulative_dist), max(xyz$cumulative_dist), col = "green2")
 # points(x, y)

  return(f(xout))
}


# Compare models for estimating centreline elevs with observations
# Creates plot
# Plot used in Bed2017 paper to support use of piece-wise linear method where observations available
# sparse_obs = 0/1 # 0 = obs have gaps | 1 = obs available for every centreline node
centreline_model_comparison<-function(path="C:/Github/synthetic_channels/meshing/development_scripts/", sparse_obs=0){

  print("****************************")
  print("AXIS LIMITS ARE HARDWIRED!!!")
  print("****************************")
  # Get test data (representative of xxxx (see Fig. in Williams et al., 2017a) where centreline elevations are extracted from parabolas 
  # constructed using known edge elevations and OMG observation values)
  elevs_fid_df=read.csv(paste0(path, "sample_centreline_points_using_OMG_for_plotting.csv"))
  seed_xyz=data.frame(x=-717233,y=-1326743,z=-158)
  mouth_xyz=data.frame(x=-700582,y=-1358375,z=-920)

  # sparse version - use to see how well methods work with gaps...
  if (sparse_obs==0){
    elevs_fid_df_SPARSE_1=elevs_fid_df[3:15,] 
    elevs_fid_df_SPARSE_2=elevs_fid_df[54:63,] 
    elevs_fid_df_SPARSE_3=elevs_fid_df[100:115,] 
    elevs_fid_df=rbind(elevs_fid_df_SPARSE_1,elevs_fid_df_SPARSE_2,elevs_fid_df_SPARSE_3)
  }

  xyz=combine_seed_mouth_cl_elevs(elevs_fid_df, seed_xyz, mouth_xyz, plotting=TRUE, verbose=verbose)

  seq(from=0, to=max(xyz$cumulative_dist), by=100)
  dists=seq(from=0, to=max(xyz$cumulative_dist), by=100)
  piece_lin_checks=matrix(NA, nrow=(length(dists)), ncol=1)
  linear_checks=matrix(NA, nrow=(length(dists)), ncol=1)
  spline_checks_fmm=matrix(NA, nrow=(length(dists)), ncol=1)
  spline_checks_periodic=matrix(NA, nrow=(length(dists)), ncol=1)
  spline_checks_natural=matrix(NA, nrow=(length(dists)), ncol=1)
  

  ########
  # linear model coeffs
  lm_out=linear_model_elevation_coeffs(elevs_fid_df, seed_xyz, mouth_xyz, plotting=TRUE) # should take in oriognal seed and mouth xyz?
  lm_coeffs=lm_out[[1]]

  for (d in 1:length(dists)){
    print(paste0("dist: ", dists[d]))

    piece_lin_checks[d]=approx(xyz$cumulative_dist, xyz$z, xout=dists[d], method="linear")$y
    linear_checks[d]=linear_model_elevations(dists[d], as.double(lm_coeffs$coefficients[1]), as.double(lm_coeffs$coefficients[2]))
    spline_checks_fmm[d]=spline(xyz$cumulative_dist, xyz$z, method = c("fmm"), xout=dists[d])$y 
    spline_checks_periodic[d]=spline(xyz$cumulative_dist, xyz$z, method = c("periodic"), xout=dists[d])$y
    spline_checks_natural[d]=spline(xyz$cumulative_dist, xyz$z, method = c("natural"), xout=dists[d])$y
  }

  piece_lin_checks<-na.omit(piece_lin_checks)
  linear_checks<-na.omit(linear_checks)

  spline_checks_fmm<-na.omit(spline_checks_fmm)
  spline_checks_periodic<-na.omit(spline_checks_periodic)
  spline_checks_natural<-na.omit(spline_checks_natural)
  
  piece_lin_df=data.frame("x"=dists, "y"=as.matrix(piece_lin_checks),"model"="Piece-wise linear") # x = cumulative dist, y = elevation
  linear_checks_df=data.frame("x"=dists, "y"=as.matrix(linear_checks),"model"="Linear") # x = cumulative dist, y = elevation
  spline_checks_fmm_df=data.frame("x"=dists, "y"=as.matrix(spline_checks_fmm), "model"="Spline (fmm)")
  spline_checks_periodic_df=data.frame("x"=dists, "y"=as.matrix(spline_checks_periodic), "model"="Spline (periodic)")
  spline_checks_natural_df=data.frame("x"=dists, "y"=as.matrix(spline_checks_natural), "model"="Spline (natural)")
  
  axis_title_size = 11
  axis_tick_size = 11
  point_label_size = 4.5
  point_line_size=1.5
  plot=ggplot(linear_checks_df,aes(x,y, colour=model)) +
    geom_line(size=point_line_size)+
    geom_line(data=spline_checks_fmm_df, aes(x,y, colour=model), size=point_line_size)+
    geom_line(data=spline_checks_periodic_df, aes(x,y, colour=model), size=point_line_size)+
    geom_line(data=spline_checks_natural_df, aes(x,y, colour=model), size=point_line_size)+
    geom_line(data=piece_lin_df,aes(x,y, colour=model), size=point_line_size) +
    geom_point(data=xyz, aes(cumulative_dist,z, colour="Elevation (obs.)"), size=2.5)+#, colour='black')+
    theme_bw() +
    xlim(c(-100,max(xyz$cumulative_dist)+100)) +
    ylim(c(-1600,0)) +
    scale_color_manual(values=c("#000000", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +
    labs(x="Distance along fjord (m)", y="Eleavation (m)") +
    theme(axis.text.x=element_text(size=axis_tick_size), 
      axis.text.y=element_text(size=axis_tick_size),
      axis.title.x=element_text(size=axis_title_size),
      axis.title.y=element_text(size=axis_title_size),
      legend.position="bottom",
      legend.text = element_text(size=axis_title_size),
      legend.title=element_blank())
    #geom_text(data=subset(xyz, cumulative_dist==min(cumulative_dist)), aes(x,y,label="Head"),size=point_label_size,check_overlap=TRUE,vjust=0.4,hjust=-0.3,nudge_x=0) +
    #geom_text(data=subset(xyz, cumulative_dist==max(min(cumulative_dist))), aes(x,y,label="Mouth"),size=point_label_size,check_overlap=TRUE,vjust=0.4,hjust=-0.3,nudge_x=0)
    
  return(plot)
}


test_approx_dealing_with_NA<-function(){

  # Create input data - contains some NA vals for z
  xyz=data.frame("z"=c(0,-100,-20,-200,-300,-500,NA,-900), "cum_dist"=c(0,100,200,300,400,500,600,700))

  # Create model output using xyz
  model=approx(xyz$cum_dist, xyz$z, method="linear")

  # Use model to calculate elevation at inout distances
  dists=seq(from=0,to=700,by=100)
  vals=approx(xyz$cum_dist, xyz$z, xout=dists, method="linear")$y

  # Plot results
  plot(dists,vals, col="red", main="Predicted (red) using model (blue)")
  lines(model$x, model$y, col="blue")

}


# Takes in the norm points with assigned elevations with use of piece-wise model
# The piecewise model will have been informed using centreline nodes (and the accumulated distance between them along the fjord)
# These centrelines can eoither be "true" - exactly halfway between the edges of a given normal, or "0" which is the mapped node 
# for which the x and y coordinates are available in out_df (cl_x,cl_y)
# This function plots elevations relative to accumulated distance using this piece-wise approach
# Exp relates to the model using centre node elevs and accum. distance
# Obs relates to the centreline elev following use o set_elevations()
#
# If "true centre" was used, the (x,y) of each centre node used to define the piece-wise will have been independent of the spacing of 
# nodes in the artificial mesh i.e. there will not be a point at eactly the "true centre" location in the point mesh - this results in Exp
# not necessarily matching Obs exactly for which the centreline node equivalent will be one of the normal nodes of as near as possible a
# distance from the mapped centre node (this was mapped using the centreline mapping algorithm). The xy of the equivalent nodes are accumulated
# and then plottied with there equivalent elevation values - any offset is related (1) to differeing distance between these nodes (than to the true 
# centre nodes) which is affected by the angle of each normal at each unique fid and (2) to the position of the node on the parabola which will 
# differ to that of the centre (at no point is the parabola flat!)
#
# If the mesh node spacing is smaller, then the equivalent node to a given true centre node will be closer and the appernt fit will be tighter.
#
# If the centreline nodes used to define the piece-wise model are equal to those that were mapped, then they have equivalent mesh nodes and thus 
# the Obs vs. Exp should match almost perfectly
#
# true_centre - if true, then the nearest node to the true centre from the mesh points is identified and used to comparse the model fit
# piecewise - if not set to 1, linear model is assumed (which results in a recalculation of the xyz$z values)
check_obs_against_centreline_model<-function(xyz, out_df, seed_xyz, mouth_xyz, true_centre=TRUE, piecewise=1, verbose=0){

  if (is.data.frame(xyz)==FALSE){xyz=as.data.frame(xyz)}
  if (is.data.frame(out_df)==FALSE){out_df=as.data.frame(out_df)}
  if (is.data.frame(seed_xyz)==FALSE){seed_xyz=as.data.frame(seed_xyz)}
  if (is.data.frame(mouth_xyz)==FALSE){mouth_xyz=as.data.frame(mouth_xyz)}

  print("Check estimated centreline elevations (will be slightly offset due to normal node spacing) match the piecewise model")

  # expected
  xyz$type="exp"

    if(piecewise!=1){
        # re-calc z using linear model - prior to this, xyz$z is the elevation of the centreline node extracted from a parabola 
        # (see combine_seed_mouth_cl_elevs() called in channel_parabola_OBS_AVAILABLE__linear_centreline())
        lm_out=linear_model_elevation_coeffs(xyz)
        lm_coeffs=lm_out[[1]]
        xyz$z=linear_model_elevations(xyz$cumulative_dist, as.double(lm_coeffs$coefficients[1]), as.double(lm_coeffs$coefficients[2]))
    }

  # observed 
  obs_df=out_df # only want the first cl_elev!!
  seed_xyz_TEMP=seed_xyz
  mouth_xyz_TEMP=mouth_xyz

  if (true_centre==TRUE){
    # get normal point for each fid that is closest to the true cl (if using the true zero approach, then this normal point is closest 
    # to the true zero and therefore should be considered when checking if the modelled elevs using the piece-wise are correct)
    #   - do this using the "dist_from_cl" field" >> if negative, it is the right hand side
    equiv_centre_x=matrix(NA, nrow=length(unique(out_df$cl_id)), ncol=1)
    equiv_centre_y=matrix(NA, nrow=length(unique(out_df$cl_id)), ncol=1)
    equiv_centre_z=matrix(NA, nrow=length(unique(out_df$cl_id)), ncol=1)
    equiv_centre_offset=matrix(NA, nrow=length(unique(out_df$cl_id)), ncol=1)

    count=0
    for (fid in  min(out_df$cl_id): max(out_df$cl_id)){

      count=count+1
      sub=subset(out_df, out_df$cl_id==fid)
      
      near_indx=which(abs(sub$cl_dist_r-abs(unique(sub$dist_from_cl)))==min(abs(sub$cl_dist_r -abs(unique(sub$dist_from_cl)))))

      sub=sub[near_indx,]

      if (unique(sub$dist_from_cl)<=0){
        # negative - right hand side
        sub=subset(sub, side==1)
      } else {
        # negative - right hand side
        sub=subset(sub, side==2)
      }

      # If more than one option, then the offset of the true to mapped centre != to the node 
      # spacing - just choose one row (the top one?) - it's only an approaximation to check 
      # the modelled elevations make sense
      equiv_centre_x[count]=sub[1,]$x
      equiv_centre_y[count]=sub[1,]$y
      equiv_centre_z[count]=sub[1,]$norm_elev
      equiv_centre_offset[count]=sub[1,]$dist_from_cl                       
    }

    equiv_centre_x<-na.omit(equiv_centre_x)
    equiv_centre_y<-na.omit(equiv_centre_y)
    equiv_centre_z<-na.omit(equiv_centre_z)
    equiv_centre_offset<-na.omit(equiv_centre_offset)

    obs_df<-data.frame("x"=matrix(equiv_centre_x),"y"=matrix(equiv_centre_y),"z"=matrix(equiv_centre_z), "equiv_centre_offset"=matrix(equiv_centre_offset))
    obs_df$type="obs"

  } else if (true_centre==FALSE) {
    
    # keep only unique rows based on col
    obs_df=subset(obs_df, !duplicated(cl_id))
    keeps=c('cl_x', 'cl_y', 'cl_elev')
    obs_df=obs_df[keeps]
    obs_df$equiv_centre_offset=0
    obs_df$type="obs"
    colnames(obs_df)=c('x','y','z','equiv_centre_offset','type')

  }
  #keep=c("true_cl_x", "true_cl_y", "cl_elev")
  #out_df_temp=out_df_temp[keep]
  #out_df_temp<-out_df_temp[!duplicated(out_df_temp),]
  #out_df_temp$type="obs"

  seed_xyz_TEMP$equiv_centre_offset=0
  seed_xyz_TEMP$type="obs"
  mouth_xyz_TEMP$equiv_centre_offset=0
  mouth_xyz_TEMP$type="obs"

  #colnames(out_df_temp)=c('x','y','z','equiv_centre_offset','type')
  #out_df_temp=rbind(seed_xyz_TEMP, out_df_temp, mouth_xyz_TEMP)
  out_df_temp=rbind(seed_xyz_TEMP, obs_df, mouth_xyz_TEMP)
  pred_df_temp=cumulative_dist(out_df_temp, verbose=verbose)

      # double check - use the cumulative values to estimate elev using the piece-wise model
      # cumulative_dist_between_centrenodes # cumulative distance between points
      # elevs_to_check                      # elevs predicted using centreline dists
    # double_check=matrix(NA, nrow=length(cumulative_dist_between_centrenodes), ncol=1)
    # for (i in 1:length(cumulative_dist_between_centrenodes)){
    #   double_check[i]=approx(xyz$cumulative_dist, xyz$z, xout=cumulative_dist_between_centrenodes[i], method="linear")$y
    # }

    # double_check_df=data.frame("cumulative_dist"=cumulative_dist_between_centrenodes, "z"=matrix(double_check), "type"="double check!")

  # Plot the model elevations (exp) with the elevations extracted from out_df (which used set_elevations(true_centre=TRUE))
  ggplot(pred_df_temp, aes(cumulative_dist, z, colour=type))+
    geom_point() +
    geom_line(data=xyz, aes(cumulative_dist, z, colour=type)) +
    #geom_point(data=double_check_df, aes(cumulative_dist, z, colour=type))+
    geom_text(data=pred_df_temp, aes(cumulative_dist, z,label=equiv_centre_offset),size=2,check_overlap=TRUE,vjust=0.4,hjust=-0.3,nudge_x=0) +
    ggtitle("Obs won't match exactly, especially if using the true centre as the elevation\nfor a given dist is relative to the normal node closest to the true\ncentre - this in turn results in a different accumulating\ndistance than considering the true centres specifically\n(not possible due to the normal node spacing)")

      
}

plot_centreline_profile<-function(norm_path, out_path, verbose=0){

  norm_df=read.csv(norm_path)
  
  # keep only unique rows based on col
  norm_df=subset(norm_df, !duplicated(cl_id))
  keeps=c('cl_x', 'cl_y', 'cl_elev')
  norm_df=norm_df[keeps]
  norm_df=cumulative_dist(norm_df, use_cl_xy=TRUE, verbose=verbose)
   
  # Plot the model elevations (exp) with the elevations extracted from out_df (which used set_elevations(true_centre=TRUE))
  cl_elev_plot=ggplot(norm_df, aes(cumulative_dist, cl_elev))+
    geom_line() +
    theme_bw() +
    xlab("Cumulative distance (m)") +
    ylab("Elevation (m)")+ 
    coord_fixed(ratio=30)

  ggsave(filename=out_path,plot=cl_elev_plot)
      
}

plot_centreline_profile_from_ras<-function(norm_path, out_path, rasF, verbose=0){

  cat("Creating plot...\n")
  norm_df=read.csv(norm_path)
  ras=raster(rasF)

  # keep only unique rows based on col
  norm_df=subset(norm_df, !duplicated(cl_id))
  keeps=c('cl_x', 'cl_y')
  norm_df=norm_df[keeps]
    
  # extract values from raster
  norm_df$smooth_elev=extract(ras, norm_df) 
  
  # calc sumulative distance between points
  norm_df=cumulative_dist(norm_df, use_cl_xy=TRUE, verbose=verbose)

  # Plot the model elevations (exp) with the elevations extracted from out_df (which used set_elevations(true_centre=TRUE))
  cl_elev_plot=ggplot(norm_df, aes(cumulative_dist, smooth_elev))+
    geom_line() +
    theme_bw() +
    xlab("Cumulative distance (m)") +
    ylab("Elevation (m)")+ 
    coord_fixed(ratio=5) #30

  cat("Saving plot...\n")
  ggsave(filename=out_path,plot=cl_elev_plot)
      
}