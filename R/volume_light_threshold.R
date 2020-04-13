# 
# library(insol)
# library(mda.lakes)
# hypso = getBathy('1881900')
# 
# lat = 43
# lon = -89.4
# 
# datetime = seq(as.POSIXct('1990-01-01'), as.POSIXct('1990-01-02'), by='min')
# 
# metjd=JD(datetime)
# 
# sunv = sunvector(metjd, lat, lon, -6)
# 
# zenith = sunpos(sunv)[,2]
# 
# Idirdif = insolation(zenith, metjd, height=300, visibility=90, RH=70,tempK=20+273.15,O3=0.02,alphag=0.1)
# 
# cos_inc_sfc=sunv%*%as.vector(normalvector(0,0)) ## or sum(sunv*normalvector(0,0))
# 
# cos_inc_sfc[cos_inc_sfc<0]=0
# Isim  = Idirdif[,1] * cos_inc_sfc


#'@title Calculate the average benthic area within light thresholds
#'
#'@param kd
#'Light attenuation value of the lake (in m^-1 units)
#'@param light_incident 
#'Vector of incident light (arbitrary units but units must match threshold units).
#'@param irr_thresh 
#'Vector of length 2 with min and max value (in that order) for light thresholds. Note: 
#'>= and <= is used for comparison, so thresholds must be adjusted accordingly.
#'@param hypso
#'Hypsography data.frame, with columns \code{depth} and \code{area}. From \code{\link{}}
#'@param area_type
#'How benthic area is to be calculated. Defaults to \code{benthic}, which tries to estimate the
#'slope corrected benthic area for each layer. \code{surface} just uses a surface area approximation, 
#'ignoring slope.
#'
#'@description
#'Returns the average amount of benthic area that is within the optical intensity thresholds based on 
#'light attenuation, hypsograph, and the desired threshold values. 
#'
#'@author Luke Winslow
#'
#'@export
area_light_threshold = function(kd, light_incident, irr_thresh=c(0,2000), hypso, area_type="benthic"){
	
  depth_area_rel <- calc_depth_area_rel(hypso, area_type)
  
	light_map <- vol_light_map(kd, light_incident, irr_thresh, hypso$depths)

	average_area <- calc_area_from_vol(light_map, depth_area_rel)
	
	return(average_area)
}


#'@title Calculate the average benthic area within temperature thresholds
#'
#'@param wtr
#'Data frame of water temperatures. See \code{\link[glmtools]{get_temp}}
#'@param wtr_thresh
#'Vector of length 2 with min and max value (in that order) for temperature thresholds. Note: 
#'>= and <= is used for comparison, so thresholds must be adjusted accordingly.
#'@inheritParams area_light_threshold
#'
#'
#'@author Luke Winslow
#'
#'@import rLakeAnalyzer
#'
#'@export
area_temp_threshold = function(wtr, wtr_thresh=c(0,25), hypso, area_type="surface"){
	
  updated_hypso <- interp_hypso_to_match_temp_profiles(wtr, hypso)
  depth_area_rel <- calc_depth_area_rel(updated_hypso, area_type)
  
  wtr_map <- vol_temp_map(wtr, wtr_thresh)
	
  average_area <- calc_area_from_vol(wtr_map, depth_area_rel)
	
	return(average_area)
}


#'@title Calculate the average benthic area within temperature and light threshold
#'
#'@inheritParams area_light_threshold
#'@inheritParams area_temp_threshold
#'
#'
#'@author Luke Winslow
#'
#'
#'@export
area_light_temp_threshold = function(wtr, kd, light_incident, irr_thresh=c(0,2000), wtr_thresh=c(0,25), hypso, area_type="surface"){
	
  updated_hypso <- interp_hypso_to_match_temp_profiles(wtr, hypso)
  depth_area_rel <- calc_depth_area_rel(updated_hypso, area_type)
  
  wtr_map <- vol_temp_map(wtr, wtr_thresh)
	light_map <- vol_light_map(kd, light_incident, irr_thresh, updated_hypso$depths)
	
	##these should theoretically be the exact same size/shape
	both_map = light_map & wtr_map  #only where both apply

	average_area <- calc_area_from_vol(both_map, depth_area_rel)
	
	return(average_area)
}

# A function that combines the separate ones in order to share code for 
#   calculating each vol map and then the areas.
area_light_temp_threshold_shared <- function(wtr, kd, light_incident, irr_thresh=c(0,2000), wtr_thresh=c(0,25), hypso, area_type="surface") {
  
  updated_hypso <- interp_hypso_to_match_temp_profiles(wtr, hypso)
  depth_area_rel <- calc_depth_area_rel(updated_hypso, area_type)
  
  light_map <- vol_light_map(kd, light_incident, irr_thresh, updated_hypso$depths)
  wtr_map <- vol_temp_map(wtr, wtr_thresh) # wtr is just daily here; upsampled below to be able to compare to light
  
  # Upsample wtr (repeat values) to match dimensions of light_map, then compare light & temp
  wtr_upsampled_map <- matrix(wtr_map, ncol = ncol(wtr_map), nrow = length(light_incident), byrow=TRUE)
  both_map <- light_map & wtr_upsampled_map  #only where both apply
  
  light_only_average_area <- calc_area_from_vol(light_map, depth_area_rel)
  temp_only_average_area <- calc_area_from_vol(wtr_map, depth_area_rel)
  light_temp_average_area <- calc_area_from_vol(both_map, depth_area_rel)
  
  return(data.frame(opti_hab = light_only_average_area, 
                    therm_hab = temp_only_average_area, 
                    opti_therm_hab = light_temp_average_area))
}

# Using the T/F map for where either light, temp, or both
#   are available, calculate the area that it equates to
calc_area_from_vol <- function(vol_map, depth_area_rel) {
  map_collapsed <- colSums(vol_map, na.rm=TRUE) # Sum number of timesteps that the condition was TRUE
  average_area <- sum(depth_area_rel * map_collapsed, na.rm=TRUE)/nrow(vol_map) # Divide by number of timesteps to get average area per day
  return(average_area)
}

vol_light_map <- function(kd, light_incident, thresholds, depths){
	
  irr_mat <- matrix(light_incident, nrow = length(light_incident), ncol = 1)
  depths_mat <- matrix(depths, nrow = 1, ncol = length(depths))
	
  light_profile <- irr_mat %*% exp(-1* kd * depths_mat) # multiple matrices together
  
  light_map = light_profile >= thresholds[1] & light_profile <= thresholds[2]
  
	#Now we need to turn it to a volumetric light map, where TRUE means 
	# that layer (not just slice) is within the thresholds
  light_map_shifted <- light_map[,-1] # Effectively "moves" values of light_map over 1 column
  light_map_to_compare <- light_map[,-ncol(light_map)] # Removes last column
  vol_light_map <- light_map_shifted & light_map_to_compare # Compares each value to the value of the next column over
	
	return(vol_light_map)
}

vol_temp_map <- function(wtr, thresholds) {
  
  wtr <- drop.datetime(wtr)
  
  # Using `drop=FALSE` to maintain matrix when there is only 1 row (otherwise becomes vector)
  wtr_map <- wtr >= thresholds[1] & wtr <= thresholds[2]
  wtr_map_shifted <- wtr_map[,-1, drop=FALSE] # Effectively "moves" values of wtr_map over 1 column
  wtr_map_to_compare <- wtr_map[,-ncol(wtr_map), drop=FALSE] # Removes last column
  vol_wtr_map <- wtr_map_shifted & wtr_map_to_compare # Compares each value to the value of the next column over
  
  return(vol_wtr_map)
}

# Moves this calculation that is used multiple times
#   into a shared function so that the code is only
#   written one time.
calc_depth_area_rel <- function(hypso, area_type) {
  if(tolower(area_type) == "surface"){
    depth_area_rel <- surface_areas(hypso$depths, hypso$areas)
  }else if(tolower(area_type) == "benthic"){
    depth_area_rel <- benthic_areas(hypso$depths, hypso$areas)
  }else{
    stop("Unrecognized area_type, must be 'surface' or 'benthic'")
  }
  return(depth_area_rel)
}

#Produces a vector of length n-1
# of benthic areas between each depth
benthic_areas <- function(depths, areas){
	
  areas_lead <- c(areas[-1], NA)
  depths_lead <- c(depths[-1], NA)
  
  trap_length <- sqrt(areas*pi) + sqrt(areas_lead*pi)
  trap_height <- sqrt( (depths_lead - depths)^2 + (sqrt(areas/pi) - sqrt(areas_lead/pi))^2 )
  benth_areas <- (trap_length * trap_height)
  benth_areas_rmna <- benth_areas[!is.na(benth_areas)]
  
  # add in area for the bottom of the lake (if comes to point, then adds 0 and doesn't change)
  last_cone <- length(benth_areas_rmna)
  benth_areas_rmna[last_cone] <- benth_areas_rmna[last_cone] + tail(areas, 1) 
  
  return(benth_areas_rmna)
}

#Produces a vector of length n-1
# 
surface_areas = function(depths, areas){
	return(-1*diff(areas))
}
