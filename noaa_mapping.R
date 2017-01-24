# Disclaimer: This is a function that will gather NOAA weather information given a dataframe of wild isolates
# and a time period to collect information for. It will calculate the average maximum, minimum, and mean during
# that time period and returns a new dataframe. Relys on input dataframe to have comments for isolation date,
# should probably be manually checked before used.


noaa_mappings <- function(df, time_period, important_trait = NULL, events = 10, 
                          additional_data = c("AA1", "AB1", "AE1", "AJ1", "AL1", "AO1", "CH1", "CI1", "CT1", "CU1", 
                                              "GA1", "GD1", "GE1", "GF1", "GH1", "GJ1", "GL1", "GM1", "GN1", "GO1", 
                                              "HL1", "IA1", "IB1", "KB1", "KE1", "KF1", "KG1", "MA1", "MF1", "MG1", 
                                              "MH1", "OB1", "RH1", "ST1")) {
  #Load required libraries
  print("Loading required libraries...")
  library(devtools)
  library(tidyr)
  library(dplyr)
  options(geonamesUsername="katiesevans")
  library(geonames)
  library(geosphere)
  library(lubridate)
  library(stationaRy)
  
  #Read wild isolate table
  wi <- df
  
  #wi with location
  wi_location <- wi %>%
    dplyr::filter(!(is.na(latitude))) %>%
    dplyr::mutate(abslat = abs(latitude), abslong = abs(longitude), elevation = c(NA))
  
  #Add elevation, NA if function fails
  print("Adding elevation...")
  for(i in 1:nrow(wi_location)) {
    wi_location$elevation[i] <- GNsrtm3(wi_location$latitude[i], wi_location$longitude[i])$srtm3
    if(wi_location$elevation[i] < 0) { wi_location$elevation[i] <- NA }
  } 
  
  #new dataframe for weather observations
  wi_weather <- wi_location %>%
    dplyr::filter(!(is.na(isolation_date))) %>%
    mutate(nearest_station = c(NA), station_distance = c(NA), start_date = c(NA), end_date = c(NA))
 
  #For monthly data, remove all strains with 'only month/year' listed
  if(time_period < 12) { 
    for(i in 1:nrow(wi_weather)) {
      if(grepl("only month/year ", wi_weather$isolation_date_comment[i], ignore.case=T)) {
        wi_weather <- filter(wi_weather, isotype != wi_weather$isotype[i])
      }
      if(grepl("only year ", wi_weather$isolation_date_comment[i], ignore.case=T)) {
        wi_weather <- filter(wi_weather, isotype != wi_weather$isotype[i])
      }
    }
    #Filter missed these somehow
    wi_weather <- filter(wi_weather, isotype != "CB4858")
    wi_weather <- filter(wi_weather, isotype != "CB4853")
    wi_weather <- filter(wi_weather, isotype != "JU397")
    wi_weather <- filter(wi_weather, isotype != "CB4857")
  } else if(time_period >= 12 && time_period < 18) {
    for(i in 1:nrow(wi_weather)) {
      if(grepl("only year ", wi_weather$isolation_date_comment[i], ignore.case=T)) {
        wi_weather <- filter(wi_weather, isotype != wi_weather$isotype[i])
      }
    }
    wi_weather <- filter(wi_weather, isotype != "CB4858")
  }
  
  #Download isd-inventory.csv to make sure data at each station is found for specific year not range
  print("Downloading NOAA station information...")
  download.file('ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-inventory.csv', 'isd-inventory.csv')
  station_inventory <- read.csv('isd-inventory.csv') %>%
    mutate(id = paste(USAF, WBAN, sep = "-"))
  
  #Download ALL weather data available
  all_data <-  additional_data
  
  #Find closest station and access data
  wi_weather_max <- wi_weather
  wi_weather_min <- wi_weather
  wi_weather_avg <- wi_weather
  #for(i in 1:nrow(wi_weather)) {
  for(i in 1:5) {
    lat = wi_weather$latitude[i]
    long = wi_weather$longitude[i]
    isolationDate <- wi_weather$isolation_date[i]
    
    #update wi_weather for start/end dates
    days <- (time_period*30)/2 #use 30 days to represent a month
    wi_weather$start_date[i] <- as.character(ymd(isolationDate) - ddays(days))
    wi_weather$end_date[i] <- as.character(ymd(isolationDate) + ddays(days))
    start <- ymd(wi_weather$start_date[i])
    end <- ymd(wi_weather$end_date[i])
    
    #Check if you need to download more than one year
    total_years <- year(end) - year(start) + 1
    st <- station_inventory %>%
      filter(YEAR >= year(start), YEAR <= year(end))
    st[st < events] <- NA #Keep years with 10+ events per month
    st <- st %>% na.omit()
    st_table <- table(st$id) %>%
      as.data.frame() %>%
      filter(Freq == total_years) #Make sure each station has correct number of years
    if(length(st_table$Var1) == 0) {
      print("Error. No station data avaiable for this time period. Entering 'NA'.")
      wi_weather$nearest_station[i] <- NA
      wi_weather$station_distance[i] <- NA
    } else {
    
      #Filter stations based on location and year with data 
      stations <- get_isd_stations(lower_lat = lat - 10,
                                   upper_lat = lat + 10,
                                   lower_lon = long - 10,
                                   upper_lon = long + 10) %>%
        mutate(dist_to_worm = c(NA)) %>%
        mutate(id = paste(usaf, wban, sep = "-")) %>%
        filter(id %in% st_table$Var1)
      
      #Find distance between all stations and worm
      for(j in 1:nrow(stations)) {
        stations$dist_to_worm[j] <- distm(c(long, lat), c(stations$lon[j], stations$lat[j]))[1]
      }
      
      #Find minimum distance station
      row <- which.min(stations$dist_to_worm)
      id <- paste(stations$usaf[row], stations$wban[row], sep="-")
      wi_weather$nearest_station[i] <- id
      wi_weather$station_distance[i] <- min(stations$dist_to_worm) / 1000 #record in km
    
      
      #Get station data using stationaRy package
      assign(paste(id, "_data", sep=""), get_isd_station_data(station_id = id, 
                                                              startyear = year(start), 
                                                              endyear = year(end),
                                                              select_additional_data = all_data))
  
      #Select all possible data (will throw an error if column not there, but thats okay)
      data_to_select <- c("wd", "ws", "ceil_hgt", "temp", "dew_point", "atmos_pres", "rh", "aa1_1", "aa1_2", "ab1_1", "ae1_1", 
                          "ae1_3", "ae1_5", "ae1_7", "aj1_1", "aj1_4", "al1_1", "al1_2", "ao1_1", "ao1_2", "ch1_2", "ch1_5", 
                          "ci1_1", "ci1_4", "ct1_1", "cu1_1", "ga1_3", "gd1_4", "ge1_3", "ge1_4", "gf1_8", "gh1_1", "gh1_4", 
                          "gh1_7", "gj1_1", "gl1_1", "gm1_1", "gm1_2", "gm1_11", "gn1_1", "gn1_2", "gn1_4", "gn1_6", "gn1_8", 
                          "gn1_10", "hl1_1", "st1_2", "ia2_1", "ia2_2", "ib1_1", "ib1_4", "ib1_7", "kb1_1", "kb1_3", "ke1_1", 
                          "ke1_3", "ke1_5", "ke1_7", "kf1_1", "kg1_1", "kg1_3", "ma1_3", "mf1_1", "mg1_1", "mh1_1", "ob1_1", 
                          "ob1_2", "ob1_5", "rh1_1", "rh1_3")
      
      #Collecting data for time period
      temp_station <- get(paste(id, "_data", sep="")) %>%
        dplyr::mutate(date = ymd(paste(year, month, day, sep="-"))) %>%
        dplyr::filter(date >= start, date <= end) %>%
        dplyr::select(one_of(data_to_select), date) %>%
        as.data.frame()
      
      #If there is no data for important_variable, choose another station
      if(!is.null(important_trait)) {
        variable <- important_trait
        while(length(which(!is.na(temp_station[[variable]]))) == 0) {
          stations_new <- stations %>%
                          filter(id != wi_weather$nearest_station[i])
          
          row <- which.min(stations_new$dist_to_worm)
          id <- paste(stations_new$usaf[row], stations_new$wban[row], sep="-")
          wi_weather$nearest_station[i] <- id
          wi_weather$station_distance[i] <- min(stations_new$dist_to_worm) / 1000
          
          #Get station data using stationaRy package
          assign(paste(id, "_data", sep=""), get_isd_station_data(station_id = id, 
                                                                  startyear = year(start), 
                                                                  endyear = year(end),
                                                                  select_additional_data = all_data))
          
          #Collecting data for time period
          temp_station <- get(paste(id, "_data", sep="")) %>%
            dplyr::mutate(date = ymd(paste(year, month, day, sep="-"))) %>%
            dplyr::filter(date >= start, date <= end) %>%
            dplyr::select(one_of(data_to_select), date) %>%
            as.data.frame()
        }
      }
    
      #Remove missing values
      if("aa1_1" %in% colnames(temp_station)) { temp_station$aa1_1[temp_station$aa1_1 == 99] <- NA }
      if("aa1_2" %in% colnames(temp_station)) { temp_station$aa1_2[temp_station$aa1_2 == 9999] <- NA }
      if("ab1_1" %in% colnames(temp_station)) { temp_station$ab1_1[temp_station$ab1_1 == 99999] <- NA }
      if("ae1_1" %in% colnames(temp_station)) { temp_station$ae1_1[temp_station$ae1_1 == 99] <- NA }
      if("ae1_3" %in% colnames(temp_station)) { temp_station$ae1_3[temp_station$ae1_3 == 99] <- NA }
      if("ae1_5" %in% colnames(temp_station)) { temp_station$ae1_5[temp_station$ae1_5 == 99] <- NA }
      if("ae1_7" %in% colnames(temp_station)) { temp_station$ae1_7[temp_station$ae1_7 == 99] <- NA }
      if("aj1_1" %in% colnames(temp_station)) { temp_station$aj1_1[temp_station$aj1_1 == 9999] <- NA }
      if("aj1_4" %in% colnames(temp_station)) { temp_station$aj1_4[temp_station$aj1_4 == 99999.9] <- NA }
      if("al1_1" %in% colnames(temp_station)) { temp_station$al1_1[temp_station$al1_1 == 99] <- NA }
      if("al1_2" %in% colnames(temp_station)) { temp_station$al1_2[temp_station$al1_2 == 999] <- NA }
      if("ao1_1" %in% colnames(temp_station)) { temp_station$ao1_1[temp_station$ao1_1 == 99] <- NA }
      if("ao1_2" %in% colnames(temp_station)) { temp_station$ao1_2[temp_station$ao1_2 == 9999] <- NA }
      if("ch1_2" %in% colnames(temp_station)) { temp_station$ch1_2[temp_station$ch1_2 == 9999] <- NA }
      if("ch1_5" %in% colnames(temp_station)) { temp_station$ch1_5[temp_station$ch1_5 == 9999] <- NA }
      if("ci1_1" %in% colnames(temp_station)) { temp_station$ci1_1[temp_station$ci1_1 == 9999] <- NA }
      if("ci1_4" %in% colnames(temp_station)) { temp_station$ci1_4[temp_station$ci1_4 == 9999] <- NA }
      if("ct1_1" %in% colnames(temp_station)) { temp_station$ct1_1[temp_station$ct1_1 == 9999] <- NA }
      if("cu1_1" %in% colnames(temp_station)) { temp_station$cu1_1[temp_station$cu1_1 == 9999] <- NA }
      if("ga1_3" %in% colnames(temp_station)) { temp_station$ga1_3[temp_station$ga1_3 == 99999] <- NA }
      if("gd1_4" %in% colnames(temp_station)) { temp_station$gd1_4[temp_station$gd1_4 == 99999] <- NA }
      if("ge1_3" %in% colnames(temp_station)) { temp_station$ge1_3[temp_station$ge1_3 == 99999] <- NA }
      if("ge1_4" %in% colnames(temp_station)) { temp_station$ge1_4[temp_station$ge1_4 == 99999] <- NA }
      if("gf1_8" %in% colnames(temp_station)) { temp_station$gf1_8[temp_station$gf1_8 == 99999] <- NA }
      if("gh1_1" %in% colnames(temp_station)) { temp_station$gh1_1[temp_station$gh1_1 == 99999] <- NA }
      if("gh1_4" %in% colnames(temp_station)) { temp_station$gh1_4[temp_station$gh1_4 == 99999] <- NA }
      if("gh1_7" %in% colnames(temp_station)) { temp_station$gh1_7[temp_station$gh1_7 == 99999] <- NA }
      if("gj1_1" %in% colnames(temp_station)) { temp_station$gj1_1[temp_station$gj1_1 == 9999] <- NA }
      if("gl1_1" %in% colnames(temp_station)) { temp_station$gl1_1[temp_station$gl1_1 == 99999] <- NA }
      if("gm1_1" %in% colnames(temp_station)) { temp_station$gm1_1[temp_station$gm1_1 == 9999] <- NA }
      if("gm1_2" %in% colnames(temp_station)) { temp_station$gm1_2[temp_station$gm1_2 == 9999] <- NA }
      if("gm1_11" %in% colnames(temp_station)) { temp_station$gm1_11[temp_station$gm1_11 == 9999] <- NA }
      if("gn1_1" %in% colnames(temp_station)) { temp_station$gn1_1[temp_station$gn1_1 == 9999] <- NA }
      if("gn1_2" %in% colnames(temp_station)) { temp_station$gn1_2[temp_station$gn1_2 == 9999] <- NA }
      if("gn1_4" %in% colnames(temp_station)) { temp_station$gn1_4[temp_station$gn1_4 == 9999] <- NA }
      if("gn1_6" %in% colnames(temp_station)) { temp_station$gn1_6[temp_station$gn1_6 == 9999] <- NA }
      if("gn1_8" %in% colnames(temp_station)) { temp_station$gn1_8[temp_station$gn1_8 == 9999] <- NA }
      if("gn1_10" %in% colnames(temp_station)) { temp_station$gn1_10[temp_station$gn1_10 == 999] <- NA }
      if("hl1_1" %in% colnames(temp_station)) { temp_station$hl1_1[temp_station$hl1_1 == 999] <- NA }
      if("st1_2" %in% colnames(temp_station)) { temp_station$st1_2[temp_station$st1_2 == 9999] <- NA }
      if("ia2_1" %in% colnames(temp_station)) { temp_station$ia2_1[temp_station$ia2_1 == 999] <- NA }
      if("ia2_2" %in% colnames(temp_station)) { temp_station$ia2_2[temp_station$ia2_2 == 9999] <- NA }
      if("ib1_1" %in% colnames(temp_station)) { temp_station$ib1_1[temp_station$ib1_1 == 9999] <- NA }
      if("ib1_4" %in% colnames(temp_station)) { temp_station$ib1_4[temp_station$ib1_4 == 9999] <- NA }
      if("ib1_7" %in% colnames(temp_station)) { temp_station$ib1_7[temp_station$ib1_7 == 9999] <- NA }
      if("kb1_1" %in% colnames(temp_station)) { temp_station$kb1_1[temp_station$kb1_1 == 999] <- NA }
      if("kb1_3" %in% colnames(temp_station)) { temp_station$kb1_3[temp_station$kb1_3 == 9999] <- NA }
      if("ke1_1" %in% colnames(temp_station)) { temp_station$ke1_1[temp_station$ke1_1 == 99] <- NA }
      if("ke1_3" %in% colnames(temp_station)) { temp_station$ke1_3[temp_station$ke1_3 == 99] <- NA }
      if("ke1_5" %in% colnames(temp_station)) { temp_station$ke1_5[temp_station$ke1_5 == 99] <- NA }
      if("ke1_7" %in% colnames(temp_station)) { temp_station$ke1_7[temp_station$ke1_7 == 99] <- NA }
      if("kf1_1" %in% colnames(temp_station)) { temp_station$kf1_1[temp_station$kf1_1 == 9999] <- NA }
      if("kg1_1" %in% colnames(temp_station)) { temp_station$kg1_1[temp_station$kg1_1 == 999] <- NA }
      if("kg1_3" %in% colnames(temp_station)) { temp_station$kg1_3[temp_station$kg1_3 == 9999] <- NA }
      if("ma1_3" %in% colnames(temp_station)) { temp_station$ma1_3[temp_station$ma1_3 == 9999.9] <- NA }
      if("mf1_1" %in% colnames(temp_station)) { temp_station$mf1_1[temp_station$mf1_1 == 99999] <- NA }
      if("mh1_1" %in% colnames(temp_station)) { temp_station$mh1_1[temp_station$mh1_1 == 99999] <- NA }
      if("mg1_1" %in% colnames(temp_station)) { temp_station$mg1_1[temp_station$mg1_1 == 99999] <- NA }
      if("ob1_1" %in% colnames(temp_station)) { temp_station$ob1_1[temp_station$ob1_1 == 999] <- NA }
      if("ob1_2" %in% colnames(temp_station)) { temp_station$ob1_2[temp_station$ob1_2 == 9999] <- NA }
      if("ob1_5" %in% colnames(temp_station)) { temp_station$ob1_5[temp_station$ob1_5 == 999] <- NA }
      if("rh1_1" %in% colnames(temp_station)) { temp_station$rh1_1[temp_station$rh1_1 == 999] <- NA }
      if("rh1_3" %in% colnames(temp_station)) { temp_station$rh1_3[temp_station$rh1_3 == 999] <- NA }
      
      #Average precipitation per hour
      if("aa1_1" %in% colnames(temp_station)) {
             temp_station <- mutate(temp_station, prcp = aa1_2 / aa1_1) %>%
               dplyr::select(-aa1_1, -aa1_2) 
      }
      
      #Average snow accumulation per hour
      if("al1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, snow_acc = al1_2 / al1_1) %>%
          dplyr::select(-al1_1, -al1_2)
      } 
      
      #Average liquid precipitation per hour
      if("ao1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, liq_prcp = ao1_2 / ao1_1) %>%
          dplyr::select(-ao1_1, -ao1_2)
      }
      
      #Average radiation per minute
      if("gm1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, glob_irr = gm1_2 / gm1_1) %>%
          mutate(uvb = gm1_11 / gm1_1) %>%
          dplyr::select(-gm1_2, -gm1_1, -gm1_11)
      }
      
      #Average solar thermal radiation per minute
      if("gn1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, up_sol_rad = gn1_2 / gn1_1) %>%
          mutate(down_therm_infr = gn1_4 / gn1_1) %>%
          mutate(up_therm_infr = gn1_6 / gn1_1) %>%
          mutate(ps_radiation = gn1_8 / gn1_1) %>%
          mutate(solar_zen_angle = gn1_10 / gn1_1) %>%
          dplyr::select(-gn1_2, -gn1_1, -gn1_4, -gn1_6, -gn1_8, -gn1_10)
      }
      
      #Ground surface observation temp per hour
      if("ia2_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, ground_surf = ia2_2 / ia2_1) %>%
          dplyr::select(-ia2_2, -ia2_1)
      }
      
      #Average dew point/wet bulb temp
      if("kg1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, dewpt_wb = kg1_3 / kg1_1) %>%
          dplyr::select(-kg1_1, -kg1_3)
      }
      
      #Average air temp
      if("kb1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, avg_air_temp = kb1_3 / kb1_1) %>%
          dplyr::select(-kb1_3, -kb1_1)
      }
      
      #Average relative humidity
      if("rh1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, avg_rh = rh1_3 / rh1_1) %>%
          dplyr::select(-rh1_1, -rh1_3)
      }
      
      #Average hourly wind
      if("ob1_1" %in% colnames(temp_station)) {
        temp_station <- mutate(temp_station, wind_max_gust = ob1_2 / ob1_1) %>%
          mutate(wind_max_dir = ob1_5 / ob1_1) %>%
          dplyr::select(-ob1_2, -ob1_5, -ob1_1)
      }

      #Average all traits and add to wi_weather
      j <- 1
      wi_weather_max$station_distance[i] <- wi_weather$station_distance[i]
      wi_weather_max$nearest_station[i] <- wi_weather$nearest_station[i]
      wi_weather_max$start_date[i] <- wi_weather$start_date[i]
      wi_weather_max$end_date[i] <- wi_weather$end_date[i]
      wi_weather_min$station_distance[i] <- wi_weather$station_distance[i]
      wi_weather_min$nearest_station[i] <- wi_weather$nearest_station[i]
      wi_weather_min$start_date[i] <- wi_weather$start_date[i]
      wi_weather_min$end_date[i] <- wi_weather$end_date[i]
      wi_weather_avg$station_distance[i] <- wi_weather$station_distance[i]
      wi_weather_avg$nearest_station[i] <- wi_weather$nearest_station[i]
      wi_weather_avg$start_date[i] <- wi_weather$start_date[i]
      wi_weather_avg$end_date[i] <- wi_weather$end_date[i]
      
      for(traits in colnames(temp_station)) {
        #Add maximum
        if(!(traits %in% colnames(wi_weather_max))){
          wi_weather_max[paste0(traits, "_max")] <- c(NA)
        }
        colnum <- which(colnames(wi_weather_max) == paste0(traits, "_max"))
        maxs <- tapply(temp_station[[j]], temp_station$date, max, na.rm = TRUE)
        is.na(maxs) <- do.call(cbind,lapply(maxs, is.infinite))
        wi_weather_max[i, colnum] <- mean(maxs, na.rm = TRUE)
        if(wi_weather_max[i, colnum] == "NaN") { wi_weather_max[i, colnum] <- NA}
        
        #Add minimum
        if(!(traits %in% colnames(wi_weather_min))){
          wi_weather_min[paste0(traits, "_min")] <- c(NA)
        }
        colnum <- which(colnames(wi_weather_min) == paste0(traits, "_min"))
        mins <- tapply(temp_station[[j]], temp_station$date, min, na.rm = TRUE)
        is.na(mins) <- do.call(cbind,lapply(mins, is.infinite))
        wi_weather_min[i, colnum] <- mean(mins, na.rm = TRUE)        
        if(wi_weather_min[i, colnum] == "NaN") { wi_weather_min[i, colnum] <- NA}
        
        #Add average
        if(!(traits %in% colnames(wi_weather_avg))){
          wi_weather_avg[paste0(traits, "_avg")] <- c(NA)
        }
        colnum <- which(colnames(wi_weather_avg) == paste0(traits, "_avg"))
        avgs <- tapply(temp_station[[j]], temp_station$date, mean, na.rm = TRUE)
        is.na(avgs) <- do.call(cbind,lapply(avgs, is.infinite))
        wi_weather_avg[i, colnum] <- mean(avgs, na.rm = TRUE)
        if(wi_weather_avg[i, colnum] == "NaN") { wi_weather_avg[i, colnum] <- NA}
        j <- j + 1
      }
    }
    #Update user on progress
    print(paste(i, '/', nrow(wi_weather), sep = ""))
  }
  
  #Join min, max, and avg; convert start_date back to dates and return
  wi_weather <- join(wi_weather_min, wi_weather_max)
  wi_weather <- join(wi_weather, wi_weather_avg) %>%
    dplyr::select(-date)
  wi_weather$start_date <- ymd(wi_weather$start_date)
  wi_weather$end_date <- ymd(wi_weather$end_date)
  return(wi_weather)
}
