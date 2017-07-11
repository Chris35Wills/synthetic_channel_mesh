#########################################################
# Program: 12c_channel_parabol_EDGE_ELEVATIONS.r
#
# Calculates depths (i.e. below zero) along each provided profile using a fixed maximum valley depth
# This can (and where possible, should) be modified to estimate parabola on observations by changing your known points from the fixed maximum valley depth as below (see python code)
# 
# This example subsets data to deal with a specific line (see 03_parabola_elevations_per_line_BULK.r for multiple line processing)
#
# IMPORTANT:
# A catch exists whereby the mid depth must always be lower than the edge elevations - currently (31st March 2016), the mid elevation is forced to a lower elev if the absolute diff 
# between the centre and the lowest edge is within a trheshold - A BETTER WAY would be to ensure the parabola meets a certain form (a value defining its order) - this is ALSO partly 
# an issue relating to the edge elevations which may change if you use the GIMP data directly (see top of script 05b*.py)
#
# @author: Chris Williams 
# @date: 31/03/16
# @version: 1.2
#########################################################

library(raster)
library(FNN)
library(tools)
library(dplyr)
#source("C:/Github/synthetic_channels/meshing/sorted_scripts/channel_parabola.r") # My R plotting function wrappers
#source("C:/Github/synthetic_channels/meshing/sorted_scripts/parabola_funcs_where_obs_available.r")
source("./channel_parabola.r") # My R plotting function wrappers
source("./parabola_funcs_where_obs_available.r")

#' Subset obs, mc and bamb points based on a distance transform (distance from land/ice where values for land/ice are 0)
#' Depending on setting of type, points will be kept depending on their distance compared to the distance of a specified point (pnt)
#' which is expected to be either the seed or mouth of a channel
#'
#' @param type = 'LT' (less than) or 'GT' (greater than)
#'      LT to be set for assigning a seed elev
#'      GT to be set for assigning a mouth elev
subset_points_by_dist<-function(obs, mc, bamb, pnt, dist, type='LT'){

  
  coordinates(mc)<- ~x+y
  coordinates(bamb)<- ~x+y
  coordinates(pnt)<- ~x+y
  crs(mc)=dist@crs
  crs(bamb)=dist@crs
  crs(pnt)=dist@crs

  mc$dist=extract(dist, mc)
  bamb$dist=extract(dist, bamb)
  pnt$dist=extract(dist, pnt)

  as.data.frame(mc)
  as.data.frame(bamb)
  as.data.frame(pnt)

  if (nrow(obs)>1){
    coordinates(obs)<- ~x+y
    crs(obs)=dist@crs
    obs$dist=extract(dist, obs)
    as.data.frame(obs)
  } 
  
  if (type=="LT"){
    print("Keeping points closer to ice/land (SEED)")
    if (nrow(obs)>1){
      obs=subset(obs, dist<=pnt$dist)
    }
    mc=subset(mc, dist<=pnt$dist)
    bamb=subset(bamb, dist<=pnt$dist)
  } else if (type=="GT"){
  
    #Only consider observations!!!! If you use the nearest bed DEM elevation or MC value, it will be out of the ocean and will give positive values!
    print("Keeping points further from ice/land (MOUTH)")
    if (nrow(obs)>1){
      obs=subset(obs, dist>=pnt$dist)
    } else {
      stop("NO OBSERVATION INFO PROVIDED -- SOME OBSERVATIONS MUST BE PROVIDED AT THE MOUTH OF THE FJORD TO TIE THE SYNTHETIC CENTRELINE ELVATIONS TO.")
    }
    #mc=subset(mc, dist>=pnt$dist)
    #bamb=subset(bamb, dist>=pnt$dist)
    mc=data.frame('x'=NA,'y'=NA,'z'=NA,'dist'=NA)
    bamb=data.frame('x'=NA,'y'=NA,'z'=NA,'dist'=NA)
  }

  return(list(obs, mc, bamb))
}

#' Find nearesst neighbour from xy2 to xy1
#' min_dist = distance of closest value in xy2 to xy1
#' min_indx = index of closest value in xy2 to xy1
#' @param details - if true returns distance and index 
#'                - if false returns index
nearest_neighbour_1stnn<-function(xy_1, xy_2){
    nn=get.knnx(xy_1, xy_2, k=1, algorithm=c("kd_tree")) # output will be a list the same length as df_obs
    min_dist=min(nn$nn.dist)       # value of shortest distance
    #min_indx=which.min(nn$nn.dist)
    min_indx=nn$nn.index[nn$nn.dist==min_dist]#which.min(nn$nn.dist)
    
    return(list(min_dist, min_indx))
}



#' Calculates the seed or mouth xyz by getting the nearest observation 
#' to the head or mouth of the fjord centreline. 
#'
#' *******
#' Could be modified to consider observations 
#' first, then mass consevtion elevs and then from Bamber2013
#' *******
#'
#' @param type - 'mouth' or 'seed' (changes which data are checked) 
#'     'mouth' - only considers observations in the ocean further from the ice sheet than the final centreline vertex
#'     'seed'  - considers any point closer to the ice sheet than the first centreline vertex
#'     'dist'  - distance transform of ocean distance from land (this should be prior calculated using the aoi land mask)
#'     'get_deepest' - if TRUE then the deepest bathymetric observation within 2km is used as opposed to just the closest
get_nearest_elev__SEED_MOUTH<-function(norm_path, dist, obs, mc, bamb, type='mouth', plotting=FALSE, get_deepest=FALSE){

  # open points
  #path_xyz=paste0(normal_filenames_path, "densified_path_0001_clipped_normals_avg_6_pnts_1200m_REARRANGED_CLIPPED.csv")
  #obs_xyz="C:/analysis_outputs_TEMP/test/vlee_79N_500m/obs_xyz.csv"
  #mc_xyz="C:/analysis_outputs_TEMP/test/vlee_79N_500m/morlighem_MC_points.csv"
  #bamber_xyz="C:/analysis_outputs_TEMP/test/vlee_79N_500m/bamber2013_gt2000m_mc.csv"

  xyz=read.csv(norm_path)   
  
  if (nrow(obs)==1){obs=data.frame('x'=NA,'y'=NA,'z'=NA)}
  if (nrow(mc)==1){mc=data.frame('x'=NA,'y'=NA,'thick'=NA,'gimp'=NA,'bed_elev'=NA)}
  if (nrow(bamb)==1){bamb=data.frame('x'=NA,'y'=NA,'z'=NA)}

  # Get point from centreline to consider
  if (type=='seed'){
      pnt=data.frame("x"=xyz[1,]$cl_x,"y"=xyz[1,]$cl_y) # point to consider: first centreline xy in the table
  } else if (type=='mouth'){
      pnt=data.frame("x"=xyz[ nrow(xyz),]$cl_x,"y"=xyz[ nrow(xyz),]$cl_y) # point to consider: last centreline xy in the table
  } 

  if (type=='seed'){
      # Consider only points closer to the land/ice
      sub_out=subset_points_by_dist(obs, mc, bamb, pnt, dist, type='LT')
      obs=sub_out[[1]]
      mc=sub_out[[2]]
      bamb=sub_out[[3]]
  } else if (type=='mouth'){
      # Consider only points away from the land/ice
      sub_out=subset_points_by_dist(obs, mc, bamb, pnt, dist, type='GT')
      obs=sub_out[[1]]
      mc=sub_out[[2]]
      bamb=sub_out[[3]]
  } else {
    stop("@param type must be either 'seed' or 'mouth'")
  } 

  # Obs
  obs_xy=data.frame('x'=obs$x, 'y'=obs$y)  
  if (nrow(obs_xy)>1 && is.na(obs$x[1])==FALSE){
    #####################
    if ( (type=='seed') | (type=='mouth' & get_deepest==FALSE)){

      print("Type is seed OR type is mouth and get_deepest is FALSE")

      nn_out=nearest_neighbour_1stnn(obs_xy, pnt) # nn_out[1] = dist nn-out[2] = index
      xyz_obs=data.frame('x'=obs[nn_out[[2]],]$x,'y'=obs[nn_out[[2]],]$y,'z'=obs[nn_out[[2]],]$z,'dist'=nn_out[[1]]) 
      
    } else if ((get_deepest==TRUE) & (type=='mouth')){

      print("Assigning deepest bathymetric observation within 2km to the mouth rather than just the closest...")
      
      nn=get.knnx(obs_xy, pnt, k=2, algorithm=c("kd_tree")) # output will be a list the same length as df_obs
      min_dists=nn$nn.dist[nn$nn.dist<2000] # value of shortest distances below a threshold dist

      if (length(min_dists)==0){
        
        print("Nearest observation > 2km from channel end - taking next nearest point...")
        nn_out=nearest_neighbour_1stnn(obs_xy, pnt) # nn_out[1] = dist nn-out[2] = index
        xyz_obs=data.frame('x'=obs[nn_out[[2]],]$x,'y'=obs[nn_out[[2]],]$y,'z'=obs[nn_out[[2]],]$z,'dist'=nn_out[[1]])         
           
      } else {
        
        min_indxs=nn$nn.index[nn$nn.dist==min_dists]

        depths=data.frame(obs[min_indxs,])
        depths$dist=min_dists
        deepest=min(depths$z)
        xyz_obs=depths[depths$z==deepest,]
        keeps=c('x','y','z', 'dist')
        xyz_obs=xyz_obs[keeps]
      }

    }

    #####################
  } else {
    xyz_obs=data.frame('x'=NA,'y'=NA,'z'=NA,'dist'=NA)
  }

  # MC
  mc_xy=data.frame('x'=mc$x, 'y'=mc$y)
  if (nrow(mc_xy)>1 && is.na(mc_xy$x[1])==FALSE){
    nn_out=nearest_neighbour_1stnn(mc_xy, pnt) # nn_out[1] = dist nn-out[2] = index
    #xyz_mc=data.frame('x'=mc[nn_out[[2]],]$x,'y'=mc[nn_out[[2]],]$y,'z'=mc[nn_out[[2]],]$bed_elev,'dist'=nn_out[[1]])
    xyz_mc=data.frame('x'=mc[nn_out[[2]],]$x,'y'=mc[nn_out[[2]],]$y,'z'=mc[nn_out[[2]],]$z,'dist'=nn_out[[1]])
  } else {xyz_mc=data.frame('x'=NA,'y'=NA,'z'=NA,'dist'=NA)}

  # Bed
  bamb_xy=data.frame('x'=bamb$x, 'y'=bamb$y)  
  if (nrow(bamb_xy)>1 && is.na(bamb_xy$x[1])==FALSE){
    nn_out=nearest_neighbour_1stnn(bamb_xy, pnt) # nn_out[1] = dist nn-out[2] = index
    xyz_bamb=data.frame('x'=bamb[nn_out[[2]],]$x,'y'=bamb[nn_out[[2]],]$y,'z'=bamb[nn_out[[2]],]$z,'dist'=nn_out[[1]])
  } else {xyz_bamb=data.frame('x'=NA,'y'=NA,'z'=NA,'dist'=NA)}

  obs_mc_bamb_xyz=na.omit(rbind(xyz_obs, xyz_mc, xyz_bamb))
  
  nearest_xyz=obs_mc_bamb_xyz[obs_mc_bamb_xyz$dist==min(obs_mc_bamb_xyz$dist),]

  if (plotting==TRUE){
    map.p <- rasterToPoints(dist)
    map.df <- data.frame(map.p)
    colnames(map.df)<-c('x','y','Distance')

    xyz_obs$type="Obs"
    xyz_mc$type="MC"
    xyz_bamb$type="Quasi"
    pnt$type=type
    nearest_xyz$type="Nearest"
    keeps=c('x','y','type')
    xyz_to_plot=na.omit(rbind(pnt[keeps], xyz_obs[keeps], xyz_mc[keeps], xyz_bamb[keeps], nearest_xyz[keeps]))
    
    map_plot = ggplot(map.df, aes(x,y)) +
      geom_raster(aes(fill=Distance)) +
      geom_point(data=pnt, aes(x,y, colour=type), pch=24, na.rm=TRUE) +
      geom_point(data=xyz_to_plot, aes(x,y, colour=type), pch=20, na.rm=TRUE) +
      theme_bw() +
      coord_equal() +
      xlab("Easting") +
      ylab("Northing")
  
    ggsave(filename=paste0(dirname(mask), "/", file_path_sans_ext(basename(norm_path)), '__nn_to_', type, '.png'),
                    plot=map_plot) 
  }
  
  keeps=c('x','y','z')
  return(nearest_xyz[keeps])

}

read_me_df<-function(data_f){
  data=read.csv(data_f)
  df<-as.data.frame(data)
  return(df)
}

check_file_exists<-function(file_to_check){
  
  existance=file.exists(file_to_check)

  if (existance==FALSE){
    stop(paste0(file_to_check, " doesn't exist!"))
  }

}

#' mask=path to raster 
#' obs_path=path to csv OR empty 
#' mc_path=path to csv 
#' bamber_path=path to csv 
#' dist_threshold = observtions within this distance are associated with a given normal
#'
#' norm_path      - the path top the normal files to use
#' norm_files     - a list of normal files to use
#' edge_path      - the path top the edge files to use - the edge files in this dir should match the normal files i.e. file_path_sans_ext(norm_file, "_edge_elevs.csv")
#' mask           - a raster object of the mask for the AOI
#' bamber_path    - path to the bamber points (or the gimp/krige points -- limit this to points within he domian, otherwise this slows everything down)
#' dist           - a distance raster (distance of ocean pixels to any other pixels)
#' obs_path       - path to the bathymetry points (limit this to points within he domian, otherwise this slows everything down)
#' dist_threshold - radius of the search window around a centreline node to search for observations
#' ignore_manual  - if using the norm_path option, this either includes or excludes paths witha  number followed by 'm' e.g. 0001 vs 0001m (manually created paths)
#' seed           - if you wish to set the seed position and elevation, pass a vector of the x-coord, y_coord and elevation e.g. seed=c(222,333,-50) otherwise leave blank and the nearest point from the obs will be selected
#' mouth          - if you wish to set the mouth position and elevation, pass a vector of the x-coord, y_coord and elevation e.g. mouth=c(222,333,-50) otherwise leave blank and the nearest point from the obs will be selected
#channel_parabola_DRIVER<-function(norm_path, edge_path, mask, mc_path, bamber_path, dist, obs_path='', true_centre=FALSE, dist_threshold=1000, ignore_manual=FALSE, BedMachinev3=FALSE, verbose=0){
channel_parabola_DRIVER<-function(norm_files, edge_path, mask, mc_path, bamber_path, dist, obs_path='', true_centre=FALSE, dist_threshold=1000, ignore_manual=FALSE, BedMachinev3=FALSE, verbose=0, seed='', mouth=''){

  check_file_exists(mask)
  check_file_exists(mc_path)
  check_file_exists(bamber_path)
  check_file_exists(dist)
  dist=raster(dist)
  if(hasValues(dist)!=TRUE){stop("Distance surface doesn't have any values...")}

#  ## Get Normal point files - output of 09_mask_clip.py << only use this if you pass the norm_path to the function rather than a list of files to use
#  norm_files=Sys.glob(paste0(norm_path, "*005*REARRANGED_CLIPPED.csv")) ## make sure zero indexed otherwise won't be sorted
#  man_norm_files=Sys.glob(paste0(norm_path, "*005*m_clipped*REARRANGED_CLIPPED.csv"))
#
#  if (ignore_manual==TRUE){
#    norm_files=setdiff(norm_files, man_norm_files)
#  }

  if (length(norm_files)==0){
    stop("No normal files found!")
  }

  norm_path=paste0(dirname(norm_files[1]), "/")
  # Catch any files that cause a fail...
  problem_files=c()
  problem_filesF=paste0(norm_path, "PROBLEMS_FILES______.csv")
  fail_count=0



  ## Get Bank elevations - see script 10_get_channel_bank_elevation.r - one file for each of the normal point files (they are specific)
  #edge_files=Sys.glob(paste0(edge_path, "*REARRANGED_CLIPPED_edge_elevs.csv"))
 
  if (BedMachinev3==FALSE){ # This first loop uses data in the old format
    # Read obs in and keep only x,y,z
    if (obs_path!=''){
      obs_xyz=read.csv(obs_path)
      #obs_xyz<-obs_xyz[c('x','y','z')]
      #keeps=c('x_bamber', 'y_bamber', 'ellips_z')
      keeps=c('x', 'y', 'layer')
      obs_xyz=obs_xyz[keeps]
      colnames(obs_xyz)<-c('x','y','z')
    } else { 
      obs_xyz=data.frame()
    }

    bamb_xyz=read.csv(bamber_path) # x,y,z

    print("Trying to read in MC...")
    mc=read.csv(mc_path, header=TRUE)
    mc_xyz=mc[c("x_bamber", "y_bamber", "thickness", "gimp_surface_elev","bed_elev")]
    colnames(mc_xyz)<-c("x", "y", "thick", "gimp", "bed_elev")
    print("MC READ IN AND COLS CHANGED...")

  } else { # This second loop uses data that were previously extracted from BedMachinev3

    obs_xyz=read.csv(obs_path, header=TRUE)
    bamb_xyz=read.csv(bamber_path, header=TRUE)
    mc_xyz=read.csv(mc_path, header=TRUE)

    header_to_keep<-c('x_bamber', 'y_bamber', 'z')  
    
    obs_xyz=obs_xyz[header_to_keep]
    bamb_xyz=bamb_xyz[header_to_keep]
    mc_xyz=mc_xyz[header_to_keep]
    
    header_to_use<-c('x','y','z')
    colnames(obs_xyz)<-header_to_use
    colnames(bamb_xyz)<-header_to_use
    colnames(mc_xyz)<-header_to_use


    ############# hardwired seeds (for NE)
    obs_xyz=rbind(obs_xyz, data.frame(x=272070, y=-801637, z=-113))
  }
  

  # open mask and calculate ocean/land distance
  msk=raster(mask)
  msk[msk!=0]=NA
  #dist=distance(msk) # ocean from land
  

  for (count in 1:length(norm_files)){
  #count=1

    current_file=norm_files[count]

    tryCatch({ 

      print(capture.output(cat("Working on ", norm_files[count])))


      # Check path length
      cut_off_clid_num=16
      df_norm=read_me_df(norm_files[count])
      query_df<-subset(df_norm, !duplicated(cl_id))

      check_not_contiguous=any(diff(unique(df_norm$cl_id))!=1)


      if (check_not_contiguous==TRUE){

        cat("XXXXXXXXXXXXX :( \n")
        cat("Skipping file - cl_id values not contiguous - the masking may have cut the channel part way along \n")
        cat("XXXXXXXXXXXXX :( \n")

      } else {

        # Check noirm and edge match - log file if not!!!
        edge_file=paste0(file_path_sans_ext(norm_files[count]), "_edge_elevs.csv")
        df_edge=read_me_df(edge_file)
        edge_match=unique(df_edge$cl_id)==unique(df_norm$cl_id) ## returns false if no match

        # Only consider paths with greater than "cut_off_clid_num" unique centreline nodes
        if ((nrow(query_df)>cut_off_clid_num) & (all(edge_match)==TRUE)){

          # Get seed and mouth xyz (path by path instance)     
          if (length(seed)==3){
            seed_xyz=data.frame(x='',y='',z='')
            seed_xyz$x=seed[1]
            seed_xyz$y=seed[2]
            seed_xyz$z=seed[3]
          } else {
            seed_xyz=get_nearest_elev__SEED_MOUTH(norm_files[count], dist, obs_xyz, mc_xyz, bamb_xyz, type="seed", get_deepest=FALSE)   # add to function --> create plot and save           
          }

          if (length(mouth)==3){
            mouth_xyz=data.frame(x='',y='',z='')
            mouth_xyz$x=mouth[1]
            mouth_xyz$y=mouth[2]
            mouth_xyz$z=mouth[3]
          } else {
            mouth_xyz=get_nearest_elev__SEED_MOUTH(norm_files[count], dist, obs_xyz, mc_xyz, bamb_xyz, type="mouth", get_deepest=TRUE) # add to function --> create plot and save            
          }       

          if(seed_xyz$z>-50){seed_xyz$z=-50}   

          print(paste0("Seed xyz: ", seed_xyz$z))
          print(paste0("Mouth xyz: ", mouth_xyz$z))

          if (seed_xyz$z < -5000 | seed_xyz$z > 200){
            stop(paste0("Spurious values set for seed: ", seed_xyz$z))
          } else if  (mouth_xyz$z < -5000 | mouth_xyz$z > 0){
            stop(paste0("Spurious values set for mouth: ", mouth_xyz$z))
          }

          if(nrow(seed_xyz)==0){ 
            stop("Seed didn't get set... check channel_parabola_DRIVER()")
          }

          if(nrow(mouth_xyz)==0){ 
            stop("Mouth value didn't get set... check channel_parabola_DRIVER()")
          }

          # Piecewise - requires seed and mouth xyz
          if (true_centre==FALSE){
            fout<-capture.output(cat(strsplit(norm_files[count], "REARRANGED*")[[1]][1], "_parabola___PIECEWISE.csv", sep=""))
            cl_fout<-capture.output(cat(dirname(norm_files[count]),"/", strsplit(basename(norm_files)[count], "*clipped")[[1]][1], "_parabola___PIECEWISE.png", sep=""))
          } else if (true_centre==TRUE){
            fout<-capture.output(cat(strsplit(norm_files[count], "REARRANGED*")[[1]][1], "_parabola___PIECEWISE_TRUE.csv", sep=""))
            cl_fout<-capture.output(cat(dirname(norm_files[count]),"/", strsplit(basename(norm_files)[count], "*clipped")[[1]][1], "_parabola___PIECEWISE_TRUE.png", sep=""))
          }
          
          #channel_parabola__piecewise(norm_files[count], edge_files[count], 

          print("About to calculate parabola (piecewise)")
          print(paste0("Norm file: ", norm_files[count]))
          print(paste0("Edge file: ", edge_file))


          ## For faster debugging! :)
          #normal_points_FILE=norm_files[count]
          #edge_elevations_FILE=edge_file
          #observations_FILE=obs_xyz
          #verbose=0
          #plotting=FALSE
          #true_centre=true_centre
          #piecewise=1
          #dist_threshold=dist_threshold
          outs=channel_parabola__piecewise(norm_files[count], edge_file, 
                                      observations_FILE=obs_xyz,
                                      fout, seed_xyz, mouth_xyz, verbose=verbose, plotting=FALSE,  
                                      true_centre=true_centre, piecewise=1, dist_threshold=dist_threshold)

          # PLOT CENTRELINE ELEVATION PROFILES
          
          plot_centreline_profile(fout, cl_fout, verbose=0)
         
          cat("Output: ", fout, "\n")

        } else if (nrow(query_df)<cut_off_clid_num) { 

          cat("Path too short so moving on to the next - need at least ", cut_off_clid_num, " unique centreline IDs", "\n")

        } else if (all(edge_match)!=TRUE) {

          cat("Files skipped: \n")
          print(paste0("Norm file: ", norm_files[count]))
          print(paste0("Edge file: ", edge_file))
          cat("This is because the edge and norm cl_id values don't match - maybe check what is going on with *_get_channel_bank_elevation.py", "\n")

        }
       
      }
    
     # end of try catch
    }, 

      error=function(e){
        cat("ERROR: Something went wrong with ", basename(norm_files[count]), "\n", conditionMessage(e), "\n")
        fail_count=fail_count+1
        problem_files[fail_count]=basename(norm_files[count])
        }) 
      

  }

  write.csv(problem_files, problem_filesF, row.names=FALSE) 

}


calc_distance_transform<-function(path, mask){
  mask=raster(mask)
  mask[mask!=0]=NA
  dist=distance(mask) # ocean from land
  ofile_dist=paste0(path, "dist_mask.tif")
  writeRaster(dist, ofile_dist, overwrite=TRUE)
  return(ofile_dist)
}


############# Example run

source("./12c_channel_parabola_EDGE_ELEVATIONS_piecewise.r")

system("mkdir -p /scratch/glacio1/cw14910/R_Raster/CHANNELS")
rasterOptions(tmpdir=paste0("scratch/glacio1/cw14910/R_Raster/CHANNELS"))

path="/scratch/glacio1/cw14910/CHANNEL_MAPPING_IMPLEMENTATION_2/N/"

norm_path=path
edge_path=path
mask=Sys.glob(paste0(path, "aoi_*_mask_500m_BedmachineV3.tif"))[1]

obs_path=paste0(path,"ibcao_bed2017bathym_COMBO_BedMachinev3__single_point.csv")
mc_path=paste0(path,"MC_points__BedMachinev3.csv")
bamber_path=paste0(path,"gimp_krig_BedMachinev3.csv")

#dist=calc_distance_transform(path, mask)
dist=paste0(path, "dist_mask.tif")

## Get Normal point files - output of 09_mask_clip.py << only use this if you pass the norm_path to the function rather than a list of files to use
norm_files=Sys.glob(paste0(norm_path, "*REARRANGED_CLIPPED.csv")) ## make sure zero indexed otherwise won't be sorted

channel_parabola_DRIVER(norm_files, edge_path, mask, mc_path, bamber_path, dist, obs_path=obs_path, true_centre=FALSE, dist_threshold=1000, ignore_manual=FALSE, BedMachinev3=TRUE)

