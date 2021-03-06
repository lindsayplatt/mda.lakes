#' @title Get Hypsometry for a given lake
#' @description
#' Returns the hypsometry profile for a lake with the given ID
#' 
#' @param site_id The character ID for the requested data
#' 
#' @return
#' Data frame with columns \code{depth} and \code{area}
#' 
#' 
#' 
#' @author 
#' Luke Winslow, Jordan Read
#'   
#'  
#' @export
getBathy	<-	function(site_id){
	
	.Deprecated('get_bathy', package='lakeattributes')
	
	return(lakeattributes::get_bathy(site_id))
	
}

#' @title Get lake surface area 
#' @description
#' Returns the surface area for a lake with given ID
#' 
#' @param site_id The character ID for the requested data (can be vector of ids)
#' 
#' @return
#'  Lake surface areas in meters^2. NA if no value available
#'  
#' @details
#' Looks for given site_ids and returns values if availalbe. All site_ids 
#' are coerced to character regardless of input type.
#' 
#' 
#' @author 
#' Jordan Read, Luke Winslow
#' 
#' @examples
#' #Multiple values
#' getArea(c('10000', '6100'))
#' 
#' #this should return NA
#' getArea('asdf')
#' 
#' @export
getArea = function(site_id){
	
	.Deprecated('get_area', package='lakeattributes')
	return(get_area(site_id))
}

#' @title Get estimated lake residence time
#' @description
#' Returns the estimated residence time for a lake with the given ID
#' 
#' @param site_id The character ID for the requested data
#' 
#' @return 
#' Estimated residence time in days
#' 
#' @details 
#' TODO
#' 
#' @references 
#' TODO Data source needed
#' 
#' @author 
#' Luke Winslow, Jordan Read
#' 
#' @examples
#' TODO
#' 
#' @export
getResidenceTime	<-	local(
	{ lookup=NULL 
		
		default.RT	<-	157.2 # days
		
		function(WBIC,default.if.null=FALSE) {
			if (is.null(lookup)) { 
				cat('Caching residence time info.\n')
				fname = system.file('supporting_files/Res.time.Diebel.csv', package=packageName())
				d	<-	read.table(fname, header=TRUE, sep=',')
				lookup <<- new.env()
				for (i in 1:nrow(d)){
					lookup[[toString(d$WBIC[i])]]	<-	d$med.RT.days[i]
				}
			}
			wbic.val = lookup[[as.character(WBIC)]]

			if (is.null(wbic.val) & default.if.null==TRUE){
			  return(default.RT)
			} else if (is.null(wbic.val) & default.if.null==FALSE){
			  return(NA)
			} else {
			  return(wbic.val)
			}
		}
	}
)

#' @title Get surrounding canopy height for a given lake
#' @description
#' Returns the surrounding canopy height for a lake with the given ID
#' 
#' @param site_id The character ID for the requested data
#' @param method Canopy height estimation method [aster or landcover]
#' @param default.if.null Default value to return if canopy height is unknown
#' 
#' @return
#'  Canopy height above lake surface level in meters
#' @details
#' TODO
#' 
#' 
#' @author 
#' Luke Winslow, Jordan Read
#' 
#' 
#' @export
getCanopy = function(site_id, default.if.null=FALSE, method="landcover"){
	
	if (tolower(method) == 'aster'){
			
		fname = system.file('supporting_files/canopyht_zonal_no_zero_num5_positivehts.csv', 
												package =packageName())
		d	= read.table(fname, header=TRUE, sep=',')
		
		vals = merge(data.frame(site_id, order=1:length(site_id)), 
								 d, by='site_id', all.x=TRUE)
		vals = vals[order(vals$order),]
		
		return(vals$MEAN_HT)
		
	}else if (tolower(method) == "landcover"){
		data(canopy)
	  id = site_id
		d = subset(canopy, `site_id` == id & source == 'nlcd')
		if(nrow(d) < 1){
		  return(NA)
		}
		return(d$canopy_m)
		
	}else{
		stop('Unidentified method ', method, ' for getCanopy [aster, landcover]')
	}
}

#' @title Get light attenuation based on long-term trend scenario for a given lake
#' @description
#' Calculate the long-term values of light attenuation for a given lake based on
#' changing scenario
#' 
#' @param site_id The character ID for the requested data
#' @param years Numeric year vector across which to calculate Kd values
#' @param year.1 The numeric year bounds of the averaging
#'  (i.e., values outside of this range will not be used)
#' @param trend The percentage of Secchi increase per year. 
#' Positive is increasing secchi (e.g., 0.94 is 0.94% increase per year)
#' @param default.if.null boolean indicating if default Kd should be used if lake has no observations
#' @return
#'  light attenuation coefficient in m^-1
#'  
#' @details
#' TODO
#' 
#' 
#' @author 
#' Luke Winslow, Jordan Read
#' 
#' 
#' 
#' @export
getScenarioKd <- function(WBIC,years,year.1=1979,year.2=2011,trend=0,default.if.null=FALSE){
  #WBIC is a string
  #years is a single numeric or vector of numerics
  #year.1 and year.2 are the bounds of the averaging (i.e., values outside of this range will not be used)
  #trend is a percentage of SECCHI increase (positive number) or decrease (negative number). E.g., 0.94 is a 0.94%/yr increase in SECCHI (decrease in Kd) . 
  #default.if.null is boolean. If default.if.null==T, a default kd will be used (centered) and the trend applied
  
  if (is.na(WBIC)){stop('WBIC cannot be NA')}
  if (is.null(WBIC)){stop('WBIC cannot be NULL')}
  
  default.kd  <-	0.6983965
  
  secchiConv  <-	1.7
  fname = system.file('supporting_files/annual_mean_secchi.txt', package=packageName())
  d	<-	read.table(fname, header=TRUE, sep='\t')
  
  useI  <-	d$WBIC==WBIC
  Kd <- vector(length=length(years))
  
  if (!any(useI) & default.if.null==T){
    year.cent = year.1  #mean(c(year.1,year.2))
    secchi.mn = secchiConv/default.kd

  } else if (!any(useI) & default.if.null==F){
    
    return(NULL)    
  } else if (any(useI) & all(is.na(d[useI,'year']))){
	  warning('No date info, using beginning of simulation as pivot point.')
	  dat.WBIC <- d[useI,]
    
    year.cent <- year.1
    secchi.mn <- mean(dat.WBIC$secchi.m.mean) # mean at pivot point!
	
  } else {
    dat.WBIC <- d[useI,]
    
    yr.i = dat.WBIC$year >= year.1 & year.2 >= dat.WBIC$year
    
    year.cent <- year.1  #mean(dat.WBIC$year[yr.i]) # the pivot point!
    secchi.mn <- mean(dat.WBIC$secchi.m.mean[yr.i]) # mean at pivot point!
  }
  
    #sech = m*x+b
    m = trend*secchi.mn*0.01
    #solve for b:
    b= secchi.mn-m*year.cent
    #convert to Kd:
    
  for (i in 1:length(years)){
    Kd[i] <- secchiConv/(m*years[i]+b)
  }
    
  return(Kd)
}

#' @title Get light attenuation coefficient for a given lake
#' @description
#' Returns the light attenuation coefficient for a lake with the given ID
#' 
#' @param site_id The character ID for the requested data
#' @param default.if.null boolean indicating if default Kd should be used if lake has no observations
#' 
#' @return
#'  Light attenuation coefficient in m^-1
#' @details
#' TODO
#' 
#' 
#' @author 
#' Luke Winslow, Jordan Read
#' 
#' @examples
#' #NA returned when no site with that ID found
#' getClarity(c('6100', '10000', 'asdf'))
#' 
#' #Default can be requested as well
#' getClarity(c('6100','asdf', '10000', 'asdf'), default.if.null=TRUE)
#' 
#' 
#' @export
getClarity = function(site_id, default.if.null=FALSE){
			
	default.kd	<-	0.6983965
	
	secchiConv	<-	1.7
	
	fname <- system.file('supporting_files/annual_mean_secchi.txt', package=packageName())
	d	<-	read.table(fname, header=TRUE, sep='\t')
	
	#lookup <<- new.env()
	#unWBIC	<-	unique(as.character(d$WBIC))
	unWBIC = site_id
	kds    = rep(NA, length(site_id))
	
	for (i in 1:length(unWBIC)){
		useI	<-	d$WBIC==unWBIC[i]
		if(any(useI)){
			secchi 	<-	d$secchi.m.mean[useI]
			attenuation.coefficient	<-	secchiConv/mean(secchi,na.rm=TRUE)
			kds[i] 	<- attenuation.coefficient
		}else{
			if(default.if.null){
				kds[i] = default.kd
			}
		}
	}
	return(kds)
}


#'@title Get surface elevation for a given lake
#'@description
#'Get the elevation of a lake with given ID
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Elevation in meters
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getElevation <- local({ lookup=NULL; 
	function(WBIC) {
	  if (is.null(lookup)) { 
	    cat('Caching elevation info\n')
	    fname <- system.file('supporting_files/WI_ManagedLakes_elevation.tsv', package=packageName())
	    d <-  read.table(fname, header=TRUE, sep='\t')
	    WBIC.names= names(d[-1]) # remove first col, it is junk
	    lookup <<- new.env()
	    for(i in 1:length(WBIC.names)){
	      names.W = WBIC.names[i] # need to remove the "X"
	      lookup[[toString(substr(x=names.W,2,stop=nchar(names.W)))]] = as.numeric(levels(d[1,names.W])[1])
	    }
	  }
	  wbic.val = lookup[[as.character(WBIC)]]
	  
	  if (is.null(wbic.val)){
	    return(NA)
	  } else {
	    return(wbic.val)
	  }
}})

#'@title Get latitude and longitude for a given lake
#'@description
#'Get the center lat/lon of a lake with given ID
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Lat/lon on the WGS84 datum as a list
#' 
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'@export
getLatLon = function(site_id){
	.Deprecated('get_latlon', package='lakeattributes')
	return(get_latlon(site_id))

}
# getLatLon <- local({ lookup=NULL; function(WBIC) {
# 	if (is.null(lookup)) { 
# 		cat('Caching lat/lon info.\n')
# 		
# 		fname <- system.file('supporting_files/WI_Lakes_WbicLatLon.tsv', package=packageName())
# 		d <- read.table(fname, header=TRUE, as.is=TRUE) 
# 		
# 		lookup <<- new.env()
# 		
# 		for(i in 1:nrow(d)){
# 			lookup[[toString(d$WBIC[i])]] = c(d$LAT[i], d$LON[i])
# 		}
# 	}
# 	lookup[[WBIC]]
# }})

#'@title Get perimeter for a given lake
#'@description
#'Get the perimeter of a lake with given ID
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Perimeter in meters
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getPerim <- local({ lookup=NULL; function(WBIC) {
	if (is.null(lookup)) { 
		cat('Caching perimeter info.\n')
		
		fname <- system.file('supporting_files/wbicAreaPerim.tsv', package=packageName())
		d <- read.table(fname, header=TRUE, as.is=TRUE) 
		
		lookup <<- new.env()
		for(i in 1:nrow(d)){
			lookup[[toString(d$WBIC[i])]] = d$PERIMETER[i]
		}
	}
	lookup[[WBIC]]
}})

#'@title Get shoreline development factor for a given lake
#'@description
#'Get the shoreline development factor (SDF) of a lake with given ID
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Return a SDF (between 1 and infinity)
#'@details
#' The ratio of a a lake's observed perimeter divided by the perimeter of a 
#' circle with the same area as the lake. Cannot be less than 1. 
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getSDF	<-	function(WBIC){
	perim	<-	getPerim(WBIC)
	area	<-	getArea(WBIC)
	circle.perim	<-	2*pi*sqrt(area/pi)
	SDF	<-	perim/circle.perim
  if (length(SDF)==0){
    return(NA)
  } else {
    return(SDF)
  }

}

#'@title Get coefficient of wind drag for a given lake or given wind sheltering coefficient
#'@description
#'Get coefficient of wind drag for a lake with a given ID or a supplied wind sheltering coefficient
#'
#'@param site_id The character ID for the requested data
#'@param Wstr The wind sheltering coefficient
#'
#'@return
#' Coefficient of wind drag 
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getCD	<-	function(site_id=NULL, Wstr=NULL, method='Hondzo'){
	
	if (is.null(site_id) & is.null(Wstr)){
		stop('either WBIC or Wstr need to be defined')
	} else if (is.null(Wstr)){
		Wstr	<-	getWstr(site_id, method=method)
	}
	
	coef_wind_drag.ref	<-	0.00140
	coef_wind_drag	<-	coef_wind_drag.ref*Wstr^0.33
	return(coef_wind_drag)
}

#'@title Get wind sheltering coefficient for a given lake
#'@description
#'Get the wind sheltering coefficient of a lake with given ID
#'
#'@param site_id The character ID for the requested data
#'@param method The desired calculation method one of c('Markfort', 'Hondzo', )
#'@return
#' The wind sheltering coefficient (between 0 and 1)
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getWstr	<-	function(site_id, method='Markfort', canopy=NULL){

	lkeArea	<-	get_area(site_id)
  
  if(is.na(lkeArea)){
    return(NA)
  }
	
	if (method=='Markfort'){
		# Markfort et al. 2010
		minWstr	<-	0.0001
		if (is.null(canopy)){
			hc	<-	max(c(getCanopy(site_id),1))
		} else {
			hc	<-	canopy
		}
		
		if(is.na(hc) | is.null(hc)){
			return(NA)
		}
		
		Xt	<-	50*hc
		
		D	<-	2*sqrt(lkeArea/pi)
		if (D<Xt){
			wind.shelter	<-	minWstr
		} else {
			wind.shelter	<-	2/pi*acos(Xt/D)-(2*Xt/(pi*D^2))*sqrt(D^2-Xt^2)
		}
	} else if (method=="Hondzo") {
		lkeArea.km2 = lkeArea*1.0e-6 # to km2
		wind.shelter= 1.0 - exp(-0.3*lkeArea.km2)
	} else {
		# Markfort et al. 2010
		minWstr	<-	0.0001
		if (is.null(canopy)){
			hc	<-	max(c(getCanopy(site_id),1))
		} else {
			hc	<-	canopy
		}
		
		Xt	<-	50*hc
		
		perim	<-	getPerim(site_id)
		shelArea	<-	perim*hc*12.5 # 25% of Markfort, as sheltering is single direction...
		shelter	<-	(lkeArea-shelArea)/lkeArea
		wind.shelter	<-	max(c(shelter,minWstr))
		if (is.null(perim)){wind.shelter<-NULL}
	}
	
	return(wind.shelter)
}

#'@title Get max depth for a given lake
#'@description
#'Get the max depth for a given lake id
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Max observed depth in meters
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getZmax <- function(site_id) {
	return(get_zmax(site_id))
}
#   if (is.null(lookup)) { 
#     cat('Caching depth info.\n')
#     fname <- system.file('supporting_files/managed_lake_info.txt', package=packageName())
#     d	<-	read.table(fname, header=TRUE, sep='\t', quote="\"")
#     
#     lookup <<- new.env()
#     ft2m  <-  0.3048
#     for(i in 1:nrow(d)){
#       mean.m <- d$max.depth.ft[i]
#       lookup[[toString(d$WBIC[i])]] = mean(ft2m*mean.m)
#     }
#   }
#   out = lookup[[WBIC]]
#   
#   #try lookup up based on bathymetry if available
#   if(is.null(out)){
#     fileN	<-	system.file(paste(c('supporting_files/Bathy/', WBIC, '.bth'), collapse=''), 
#                          package=packageName())
#     if (file.exists(fileN)){
#       data	<-	read.table(fileN,header=TRUE,sep='\t')
#       out = max(data$depth, na.rm=TRUE)
#     }
#   }
#   
#   return(out)
# }})

#'@title Get mean depth for a given lake
#'@description
#'Get the mean depth for a given lake id
#'
#'@param site_id The character ID for the requested data
#'
#'@return
#' Mean calculated depth in meters
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getZmean	<-	function(WBIC){
	ft2m	<-	0.3048
	fname <-  system.file('supporting_files/managed_lake_info.txt', package=packageName())
	data	<-	read.table(fname, header=TRUE, sep='\t', quote="\"")
	useI	<-  data$WBIC==as.numeric(WBIC)
	mean.depth	<-	NA
	if (any(useI)){
		mean.depth	<-	ft2m*mean(data$mean.depth.ft[useI],na.rm=TRUE)
	}
	if (is.na(mean.depth)){
		mean.depth	<-	1/3*getZmax(WBIC)
	}
	return(mean.depth)
}

#'@title Get the ice on date for a given lake
#'@description
#'Get the ice on date for a given lake id and year
#'
#'@param site_id The character ID for the requested data
#'@param year The year for which you need the ice-on date
#'@return
#' Max observed depth in meters
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getIceOn	<-	local({ lookup=NULL; fcache=''; function(WBIC, year, fname='empirical.ice.tsv'){
  early.freeze	<-	8	# month of year
	# the ice on for each lake for a given year
	# ice on is assumed to either happen during the same calendar year, or within the NEXT year
  if (is.null(lookup) || fname != fcache) { 
  	cat('Caching ice info.\n')
  	#cache filename so we know which we've got cached
  	fcache <<- fname
  	
	  fname <- system.file(paste0('supporting_files/', fname), package=packageName())
		empir.ice = read.table(fname, sep='\t', header=TRUE, as.is=TRUE) 
	  empir.ice$posix = as.POSIXct(empir.ice$DATE)
		lookup <<- new.env()
		lookup[["dt"]] = empir.ice
  }
  
  empir.ice = lookup[["dt"]]
  
  first.plausible = as.POSIXct(paste0(year,'-08-01'))
  last.plausible = as.POSIXct(paste0(year+1,'-02-15'))
	
  ice.on	<-	vector(length=length(WBIC))
  
	for (j in 1:length(WBIC)){
		
		use.i	<- WBIC[j]==empir.ice$WBIC & 
											empir.ice$ON.OFF=="on" & 
											empir.ice$posix > first.plausible & 
											empir.ice$posix < last.plausible
		
		pos.results = empir.ice$DATE[use.i]
		
		if(length(pos.results) > 1){
			warning(sprintf("Ambiguous ice on in year %i for lake %s, using last value", year, WBIC))
			ice.on[j]<- tail(pos.results,1)
		}else{			
			ice.on[j]<- pos.results
		}
		
	}
	return(ice.on)
}})

#'@title Get the ice off date for a given lake
#'@description
#'Get the ice off date for a given lake id and year
#'
#'@param site_id The character ID for the requested data
#'@param year The year for which you need the ice-off date
#'@return A vector of ice off dates (as character %Y-%m-%d) with the same length as supplied \code{site_id}
#'
#'@details
#'TODO
#'
#'
#'@author 
#'Luke Winslow, Jordan Read
#'
#'
#'
#'@export
getIceOff	<-	function(site_id, year, fname='empirical.ice.tsv') {
	# the ice off for each lake for a given year
	# ice off is assumed to happen during the same calendar year

	fname <- system.file(paste0('supporting_files/', fname), package=packageName())
	empir.ice = read.table(fname, sep='\t', header=TRUE, as.is=TRUE) 
	
	ice.off	<-	vector(length=length(site_id))
	for (j in 1:length(site_id)){
		use.i	<-	site_id[j]==empir.ice$WBIC & empir.ice$ON.OFF=="off" & substr(empir.ice$DATE,1,4)==as.character(year)
		if (any(use.i)){
			pos.results	<-	empir.ice$DATE[use.i]
			ice.off[j]	<-	pos.results[1] # warn if length(pos.results)==2?
		} else { 
			ice.off[j]	<-	NA
		}
		
	}
	return(ice.off)
}
