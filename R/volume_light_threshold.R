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
  
	light_map = vol_light_map(kd, light_incident, irr_thresh, hypso$depths)

	light_map_collapsed = apply(light_map, 2, sum, na.rm=TRUE)

	average_area = sum(depth_area_rel * light_map_collapsed)
	
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
  
  wtr = drop.datetime(wtr)
	
	map = wtr >= wtr_thresh[1] & wtr <= wtr_thresh[2]
	
	vol_map = map[,-1] #just needs to be 1 less column than light_map
	
	for( i in 1:(ncol(map)-1) ){
		vol_map[,i] = map[,i] & map[,i+1]
	}
	
	map_collapsed = apply(vol_map, 2, sum, na.rm=TRUE)
	
	average_area = sum(depth_area_rel * map_collapsed, na.rm=TRUE)
	
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
  
	wtr = drop.datetime(wtr)
	
	map = wtr >= wtr_thresh[1] & wtr <= wtr_thresh[2]
	vol_map = map[,-1] #just needs to be 1 less column than light_map
	
	for( i in 1:(ncol(map)-1) ){
		vol_map[,i] = map[,i] & map[,i+1]
	}
	
	light_map = vol_light_map(kd, light_incident, irr_thresh, updated_hypso$depths)
	
	##these should theoretically be the exact same size/shape
	both_map = light_map & vol_map  #only where both apply

	map_collapsed = apply(both_map, 2, sum, na.rm=TRUE)
	
	average_area = sum(depth_area_rel * map_collapsed, na.rm=TRUE)
	
	return(average_area)
}


vol_light_map = function(kd, light_incident, thresholds, depths){
	
	light_profile = data.frame(Io=light_incident)
	
	for(i in 1:length(depths)){
		colname = paste0('irr_', depths[i])
		light_profile[,colname] = light_profile$Io * exp(-1* kd * depths[i])
	}
	
	#drop incident light
	light_profile = light_profile[,-1]
	
	#Now we need to turn it to a volumetric light map, where TRUE means 
	# that layer (not just slice) is within the thresholds
	light_map = light_profile >= thresholds[1] & light_profile <= thresholds[2]
	
	vol_light_map = light_map[,-1] #just needs to be 1 less column than light_map
	
	for( i in 1:(ncol(light_map)-1) ){
		vol_light_map[,i] = light_map[,i] & light_map[,i+1]
	}
	
	return(vol_light_map)
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
  
  return(benth_areas_rmna)
}

#Produces a vector of length n-1
# 
surface_areas = function(depths, areas){
	return(-1*diff(areas))
}
