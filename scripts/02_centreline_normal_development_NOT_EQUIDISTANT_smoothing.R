##############################
# Program: 02_centreline_normal_development_NOT_EQUIDISTANT_smoothing.r
#
# Same as 01_centreline_normal_development_NOT_EQUIDISTANT.r but calculates the vector to which normals 
# are relative to, based on neighbours that are further from the point within the iteration loop (i.e. 
# for point x, calculate vector using points x-n and x+n as opposed to x+1 and x-1). Not using this smoothing 
# approach makes the normals very sensitive to the fairly coarse centreline construction, making them cross 
# over each other all over the palce... >> this filtering step is very similar to that used in Goff and Nordfjord, 2004.
#
# Calculates channel width relative to the normal of a centreline point consideiring 
# the centreline preceding and succeding neighbours. Points are constructed along the 
# normal for a width of xx m. The nearest neighbour to the points along the normal 
# selected from the edge points is selected and its distance is calculated relative to 
# the centreline point. This distance is doubled to give channel width. The relationship 
# between width and length along the channel is also calculated.
#
# A csv containing xy coords of the centreline points, 2xdistance to the nn 
# (channel width) and the coords of the nn for plotting and considering if the 
# relationship between the 2 points is reasonable.
# 
# This is very helpful (considering my mathsiness...): http://stackoverflow.com/questions/7469959/given-2-points-how-do-i-draw-a-line-at-a-right-angle-to-the-line-formed-by-the-t/7470098#7470098
#
# @ Chris Williams (+help from Dom Huelse and Alba Martin) 10/11/15
##############################

library(sp)
library(pgirmess)
library(rgdal)
library(maptools)
library(tools)


vector_from_norms<-function(in_path, out_path, file_glob="densified*.csv", verbose=TRUE, dist_between_norm_vertices=100, max_length_of_side=8000){

	filenames=Sys.glob(paste0(in_path,file_glob))
	
	file_counter<-0
	for (file in filenames){

		file_counter<-file_counter+1

		file_no_ext<-file_path_sans_ext(basename(file))

		print(paste0("Working on: ", file_no_ext))

		# Read the .csv file
		data <- read.csv(file, stringsAsFactors = FALSE)
		#print(typeof(data))

		# look at the data structure
		#str(data)

		# convert list to matrix
		#mat_x<-as.matrix(data$coord_x)
		#mat_y<-as.matrix(data$coord_y)
		mat_x<-as.matrix(data[,1])
		mat_y<-as.matrix(data[,2])

		mat_xy<-cbind(mat_x, mat_y)
		remove(mat_x, mat_y) # clean up variables

		# remove duplicates
		cl_pnts_xy<-unique(mat_xy)

		#1 Get the points either side of each centre line point (starting from 0+1 and ending at end-1) <-- NB/ points should already be ordered
		#2 Calculate the vector between these points
		#3 Calculate the left and right normal vectors to that calc. in #2
		#4 Normalise #3 (to get a vector unit of 1 - you can check this remember)
		#5 Using similar (ideally the same) spacing as for the centreline points (100m?), use this scaler to calculate point locations along the normal in both directions and get the xy values
		#6 Repeat to create multiple points
		#7 Write the normal point xy values in an xy matrix - add a second set of xy columns giving the xy of the centreline point and another column with the centreline ID

		#~~~~~~~~~~~~~~~~~~~~~~~
		#spacing=100. # 100
		#max=8000 # 5000.
		spacing=dist_between_norm_vertices
		max=max_length_of_side 

		#check max%%spacing==0 otherwise fail...
		num_points=(max/spacing)
		nrow_guesstimate=1000000 # << something big!!!!

		#predefine length of lists to append to (faster):										<<<<<<<<<< do this out of the loop so you'll get one big matrix containing points for all centreline locations (can filter matrix based on FID column)
		#http://stackoverflow.com/questions/20730537/add-new-row-to-matrix-one-by-one
		#r_right_from_m_list <- matrix(NA, nrow=nrow_guesstimate, ncol=5) 
		#r_left_from_m_list <- matrix(NA, nrow=nrow_guesstimate, ncol=5)

		#r_normal_pnts_from_m_list <- matrix(NA, nrow=nrow_guesstimate, ncol=11) 
		normal_pnts_from_m_list <- matrix(NA, nrow=nrow_guesstimate, ncol=7) 

		# loop thorugh all centre points and calculate normals
		#~~~~~~~~~~~~~~~~~~~~~~~
		vector_pnt_spacing<-6 #e.g. a value of 2 will use points 2 indices before and after from which to calculate the vector
		vector_smoothing_window_m_equiv<-(vector_pnt_spacing*spacing)*2

		#strt<-(vector_pnt_spacing+1)
		#end<-(nrow(cl_pnts_xy))-(vector_pnt_spacing)
		strt<-vector_pnt_spacing+1
		end<-nrow(cl_pnts_xy)-vector_pnt_spacing


		# skip path if too short -- path is too short where the number of points followint the starting fid are less than the vector spacing - in this case, the path is skipped
		#if (end>strt){
		if (end>vector_pnt_spacing){

			for (i in strt:end) { 
			# starts at vector_pnt_spacing+1 to ensure the first point considered can be compared the points from the start of the string of points according to the value of vector_pnt_spacing
			# ends at (length(cl_pnts_xy)/2)-vector_pnt_spacing to ensure the last point considered can compare the points to the end of the string of points according to the value of vector_pnt_spacing
				
				if (verbose==TRUE){ 
					cat("Working on centre point ", i, "\n") 
					cat("start: ", strt, "\n")
					cat("end: ", end, "\n")
					}

				###### Step 1
				# Define points
				# Ensure only unique rows are used when considering the input data set
				# NB/ pnt_a is downstream of pnt_b as the largest centreline point FID number is at the start of the channel

				pnt_a = cl_pnts_xy[i-vector_pnt_spacing,] # point before
				pnt_b = cl_pnts_xy[i+vector_pnt_spacing,] # point after
				pnt_twixt = cl_pnts_xy[i,]   # point between
								
				###### Step 2
				# Calculate the directing vectors

				directing_vector_i_a_twxt = pnt_twixt[1]-pnt_a[1] # a -> twixt
				directing_vector_j_a_twxt = pnt_twixt[2]-pnt_a[2] # a -> twixt

				directing_vector_i_b_twxt = pnt_b[1]-pnt_twixt[1] # b -> twixt
				directing_vector_j_b_twxt = pnt_b[2]-pnt_twixt[2] # b -> twixt
					
				a_twxt=c(directing_vector_i_a_twxt, directing_vector_j_a_twxt)
				b_twxt=c(directing_vector_i_b_twxt, directing_vector_j_b_twxt)

				###### Step 3
				# Calculate unit vectors

				r_a_twxt_len=sqrt((directing_vector_i_a_twxt)^2+(directing_vector_j_a_twxt)^2)
				r_b_twxt_len=sqrt((directing_vector_i_b_twxt)^2+(directing_vector_j_b_twxt)^2)

				r_a_twxt_i_unit=directing_vector_i_a_twxt/r_a_twxt_len
				r_a_twxt_j_unit=directing_vector_j_a_twxt/r_a_twxt_len

				r_b_twxt_i_unit=directing_vector_i_b_twxt/r_b_twxt_len
				r_b_twxt_j_unit=directing_vector_j_b_twxt/r_b_twxt_len
					
				###### Step 4
				# Add unit vectors to create a summed vector "w"

				w_i=r_a_twxt_i_unit+r_b_twxt_i_unit
				w_j=r_a_twxt_j_unit+r_b_twxt_j_unit

				###### Step 5
				# Calculate the Left and Right normals of w
					
				left_norm_x = (w_j*-1)
				left_norm_y =  w_i
				left_norm = c(left_norm_x, left_norm_y)

				right_norm_x =  w_j
				right_norm_y = (w_i*-1)
				right_norm = c(right_norm_x, right_norm_y)

				###### Step 6
				# Calculate unit vectors of L and R normals

				w_left_len=sqrt((left_norm_x)^2+(left_norm_y)^2)
				w_right_len=sqrt((right_norm_x)^2+(right_norm_y)^2)

				w_left_i_unit=(1/w_left_len)*left_norm_x
				w_left_j_unit=(1/w_left_len)*left_norm_y

				w_right_i_unit=(1/w_right_len)*right_norm_x
				w_right_j_unit=(1/w_right_len)*right_norm_y

				###### Step 7
				##calc a single point along both normals
				#spacing=0.5

				#pnt_norm_right_x=(w_right_i_unit*spacing)+pnt_twixt[1]
				#pnt_norm_right_y=(w_right_j_unit*spacing)+pnt_twixt[2]
				#pnt_norm_right_xy = c(pnt_norm_right_x, pnt_norm_right_y)
				##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#pnt_norm_left_x=(w_left_i_unit*spacing)+pnt_twixt[1]
				#pnt_norm_left_y=(w_left_j_unit*spacing)+pnt_twixt[2]
				#pnt_norm_left_xy = c(pnt_norm_left_x, pnt_norm_left_y)


				###### Step 7
				#create 2km worth of points either site of pnt_twixt
				#create points every 100m
				
				#for (ii in 1:num_points){ # <<< this overwrites the first 1->num_points rows of the left and right point matrixes each time
				
				i_strt=i-vector_pnt_spacing # this assumes i starts at 2 and is only used to ensure the top columns of the outoput list are populated (if you don't believe me, replace "i_strt" with i and see what happens)
				#spacing=100.

				for (ii in (((i_strt-1)*num_points)+1):(((i_strt-1)*num_points)+num_points)){ # gives the row inidces of the output matrix - jiggery pokery is to make sure you loop through the correctly alloted lines
					
					########################
					## Calculate steps from centreline (independent of output index)
					strt<-(((i_strt-1)*num_points)+1)-1
					max_steps<-(((i_strt-1)*num_points)+num_points) # index number of max steps in this loop cycle
					steps_left<-(ii-max_steps)*-1
					steps_total<-max_steps-strt
					steps_made<-steps_total-steps_left
			 		max_steps<-(((i_strt-1)*num_points)+num_points) 
					########################

					dist=steps_made*spacing # this should range by the number of steps and NOT the row indices
					#print(dist)

					#r_right_i=(r_right_norm_i*dist)+pnt_twixt[2]
					#r_right_j=(r_right_norm_j*dist)+pnt_twixt[1]

					#r_left_i=(r_left_norm_i*dist)+pnt_twixt[2]
					#r_left_j=(r_left_norm_j*dist)+pnt_twixt[1]

					pnt_norm_right_x=(w_right_i_unit*dist)+pnt_twixt[1]
					pnt_norm_right_y=(w_right_j_unit*dist)+pnt_twixt[2]
					
					pnt_norm_left_x=(w_left_i_unit*dist)+pnt_twixt[1]
					pnt_norm_left_y=(w_left_j_unit*dist)+pnt_twixt[2]

					#append to lists...
					#r_normal_pnts_from_m_list[ii,]<-c(pnt_norm_right_x, pnt_norm_right_y, pnt_norm_left_x, pnt_norm_left_y,pnt_twixt[1],pnt_twixt[2],i,pnt_a[1],pnt_a[2],pnt_b[1],pnt_b[2]) # nb/ i in the 7th colum is the position in the array of the centreline points NOT necessarily the FID
					normal_pnts_from_m_list[ii,]<-c(pnt_norm_right_x, pnt_norm_right_y, pnt_norm_left_x, pnt_norm_left_y,pnt_twixt[1],pnt_twixt[2],i) # nb/ i in the 7th colum is the position in the array of the centreline points NOT necessarily the FID
				
				}
			}

			# Remove NaN values from dataframe
			#colnames(r_normal_pnts_from_m_list)<-c("x_right","y_right","x_left","y_left","cl_x","cl_y","cl_id","a_x","a_y","b_x","b_y")
			normal_pnts_from_m_list<-na.omit(normal_pnts_from_m_list) # remove NaNs

			##################
			## Calculate distance of each point from its centre
			cat("Calculating distances of points to centre...\n")
			
			euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

			cl_xy<-normal_pnts_from_m_list[,5:6]

			# right norm
			r_xy<-normal_pnts_from_m_list[,1:2]
			dist_r <- NULL
			for(i in 1:nrow(normal_pnts_from_m_list)) dist_r[i] <- euc.dist(r_xy[i,],cl_xy[i,])
			dist_r

			# left norm
			l_xy<-normal_pnts_from_m_list[,3:4]
			dist_l <- NULL
			for(i in 1:nrow(normal_pnts_from_m_list)) dist_l[i] <- euc.dist(l_xy[i,],cl_xy[i,])
			dist_l

			path<-file_counter

			#normal_pnts_from_m_list<-cbind(normal_pnts_from_m_list, dist_r, dist_l)
			normal_pnts_from_m_list<-cbind(normal_pnts_from_m_list, dist_r, dist_l, path)
			#colnames(normal_pnts_from_m_list)<-c("x_right","y_right","x_left","y_left","cl_x","cl_y","cl_id","cl_dist_r","cl_dist_l")
			colnames(normal_pnts_from_m_list)<-c("x_right","y_right","x_left","y_left","cl_x","cl_y","cl_id","cl_dist_r","cl_dist_l","path")

			fout=capture.output(cat(out_path,"/", file_no_ext, "_normals_avg_", vector_pnt_spacing,"_pnts_", vector_smoothing_window_m_equiv, "m.csv", sep=""))
			#fout=capture.output(cat("C:/Github/synthetic_channels/test_output/test_meshing/r_left_and_right_from_m_list_non_equidistant_vector_avg_", vector_pnt_spacing,"_pnts_", vector_smoothing_window_m_equiv, "m_TEST_KINKED.csv", sep=""))
			dir.create(dirname(fout), showWarnings=FALSE) # create folder if it doesn't exist

			cat("Writing to file...\n") 
			write.csv(normal_pnts_from_m_list, file=fout, row.names=FALSE)
		}
	}

	cat("****************\n")
	cat("Meshing complete.\n")
	cat("****************\n")
}


in_path="../test_outputs/"
out_path=in_path
file_glob="densified_*.csv"
vector_from_norms(in_path, out_path, file_glob=file_glob, dist_between_norm_vertices=50, verbose=FALSE)