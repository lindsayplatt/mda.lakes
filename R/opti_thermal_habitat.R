
#'@title Calculate optical thermal habitat using temp and light thresholds 
#'
#'@inheritParams area_light_temp_threshold
#'@param interp_hour Interpolate to hourly data (and )
#'
#'@return data.frame with three columns. opti_hab, therm_hab, opti_therm_hab 
#'for areas of each habitat type (with opti_therm_hab being the overlap of both). Units are in
#'m^2*days. Divide by the number of days in the data to 
#'
#'
#'@export
opti_thermal_habitat = function(opt_wtr, io, kd, lat, lon, hypso, irr_thresh=c(0,2000), wtr_thresh=c(0,25), area_type="benthic", interp_hour=FALSE){
	
# 	nml = read_nml(nml_file)
# 	
# 	kd = get_nml_value(nml, 'Kw')
# 	lat = get_nml_value(nml, 'latitude')
# 	lon = get_nml_value(nml, 'longitude')
# 	
# 	#get hypso data and interp to standard resolution
# 	hypso = get_hypsography(nml)
# 	names(hypso) = c('depth', 'area')
# 	hypso = interp_hypso(hypso, dz=0.1, force_zero_area = TRUE)
# 	
	#also extract temp at the same resolution as hypsography. This will help us standardize across light and temp
	wtr = opt_wtr
	
#	io = get_var(nc_file, 'I_0')
	
	#Now, if we are going to interp this, we need to interp io and wtr to at least hourly or so
	if(interp_hour){
	  io <- create_irr_day_cycle(lat,lon, dates=io[[1]], irr_mean = io[[2]], by='min')
	  habitat_shared <- area_light_temp_threshold_shared(wtr, kd, io[[2]], irr_thresh, wtr_thresh, hypso, area_type)
	}else{
		light_alone = area_light_threshold(kd, io[,2], irr_thresh, hypso, area_type)
		temp_alone  = area_temp_threshold(wtr, wtr_thresh, hypso, area_type)
		light_temp  = area_light_temp_threshold(wtr, kd, io[,2], irr_thresh, wtr_thresh, hypso, area_type)
		habitat_shared <- data.frame(opti_hab=light_alone, therm_hab=temp_alone, opti_therm_hab=light_temp)
	}
	
	
	
	return(habitat_shared)
}


#'@title Linear interpolate hypsometry 
#'
#'@description Interpolate bathymetry to a specified depth interval. Includes 
#'options to ensure validity of profile (force zero area layer, monotonically 
#'decreasing area)
#'
#'@param hypso Hypsography data.frame with columns \code{depth} (meters, positive downward) and \code{area} (m^2)
#'@param dz Requested resample depth interval
#'@param force_zero_area Force the hypsography profile to keep a zero-area depth, even if the input did not have
#'a zero-area depth, or if the zero area depth does not fall perfectly on the \code{dz} interval
#'
#'@return A hypso data.frame with columns \code{depth and area}
#'
#'@examples
#'hypso = getBathy('439800')
#'new_hypso = interp_hypso(hypso, dz=0.1)
#'
#'plot(hypso$area, hypso$depth, ylim=rev(range(new_hypso$depth)), type='o')
#'lines(new_hypso$area, new_hypso$depth, lty=2, col='red')
#'
#'@export
interp_hypso = function(hypso, dz=0.1, force_zero_area=TRUE){
	
	#check data frame has what we need
	if(!all(names(hypso) %in% c('areas', 'depths'))){
		stop('Bathy data.frame must have columns areas and depths (currently: ', paste(names(bathy), collapse=","), ')')
	}
	
	new_hypso = data.frame( depths=seq(min(hypso$depths), max(hypso$depths), by=dz) )
	
	
	new_hypso$areas  = approx(hypso$depths, hypso$areas, new_hypso$depths)$y
	
	##Now we have four cases
	#1) We don't want to force zero hypso area
	#2) OR new_hypso happens to have a zero value. Great, return!
  #return whatever we have
	if(!force_zero_area || any(new_hypso$areas == 0)){
		return(new_hypso)
	}
	
	
	#3) new_hypso lost the zero that was on the hypso supplied
	if(!any(new_hypso$areas == 0) && any(hypso$areas == 0) ){
		new_hypso = rbind(new_hypso, hypso[hypso$areas==0,])
		new_hypso = new_hypso[order(new_hypso$depths), ]
	}
	
	#3) new_hypso lost the zero that was on the hypso supplied
	# in this case, add dz to max depth, and call that zero
	# Maybe TODO: Extrapolate a zero downward linearly?
	if(!any(new_hypso$areas == 0) && !any(hypso$areas == 0) ){
		new_hypso[nrow(new_hypso)+1, 'depths'] = max(new_hypso$depths) + dz
		new_hypso[nrow(new_hypso), 'areas'] = 0
		new_hypso = new_hypso[order(new_hypso$depths), ]
	}
	return(new_hypso)
	
}

interp_hypso_to_match_temp_profiles <- function(wtr, hypso) {
  
  # Match hypso depths to water temperature profile depths
  matched_depths <- get.offsets(wtr)
  matched_areas <- approx(hypso$depths, hypso$areas, xout=matched_depths)$y
  matched_hypso <- list(depths = matched_depths, areas = matched_areas)
  
  return(matched_hypso)
}
