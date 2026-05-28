#==== daylength data and fit to 3 parat Anatone =====
dl <- read.csv("data/anatone.daylengths.csv")
daylength <- cbind.data.frame(day=1:365,time=NA)
doy <- 0
for(i in 2:13){for(j in 1:31){
  if(!is.na(dl[j,i])){
    doy <- doy + 1; 
    daylength$day[doy] <- doy;
    hold <- as.numeric(str_split(dl[j,i],":")[[1]]); 
    daylength$time[doy]<- hold[1] + hold[2]/60
  }
}}


#==== Water Flow & Temp. data ====
# Flow itself is not used in the current model

water0 <- read.csv("./data/usgs.anatone.data.csv")
head(water0)

#==== Air data ====

# Lewiston Air 
# Request was made at web site behind https://www.ncdc.noaa.gov/cdo-web/
# site is Lewiston,ID
# data set is 'Daily Summaries'
# data through a complete water year (ends 9-30)

# ONLY AIR and PRECIP 
air0 <- read.csv("./data/Lewiston.Air.csv")

air0 <- as_tibble(air0) %>% mutate(Date=as.Date(DATE,format="%m/%d/%Y")) %>% 
  mutate(TavgComp = (TMAX + TMIN)/2 ) %>% mutate(Precip=PRCP) %>%
  dplyr::select(-c(STATION,NAME,DATE,TAVG,TMAX,TMIN,SNOW,SNWD)) 
air0$TavgComp <- zFiller(air0$TavgComp)
air0 <- air0[!is.na(air0$TavgComp),]
air0 <- air0 %>% mutate(year = year(Date)) %>% mutate(jul=yday(Date))%>% mutate(Tair = (TavgComp-32)*5/9  ) 
# Move everything to  Water Year time. Clip for 365 day year and end WY 2022

air1 <- cbind.data.frame(air0,zDOY.WY(cbind.data.frame(air0$year,air0$jul))[,3:4])
air1 <- air1[air1$WY >= 1961 & air1$WY <= 2022,]
air1 <- air1[air1$WYdoy <= 365,]
air1$Precip[is.na(air1$Precip)] <- 0 

air2 <- air1  %>% rename(doy=WYdoy) %>% 
  dplyr::select(-year) %>% 
  rename(year=WY) %>% 
  dplyr::select(-PRCP) %>%
  dplyr::select(-Precip)

df <- left_join(water0,air2 %>% dplyr::select(-Date),by=c("year","doy")) 

#==== Air regime metrics  ====

metrics <- NULL
for(y in unique(air2$year)){
  I <- air2$year == y
  temps <- df$ana.temp[df$year==y]
  air <- air2$Tair[I]
  springairA <-  mean(air[183:212],na.rm=TRUE)
  springairMar <-  mean(air[152:182],na.rm=TRUE)
  springairMay <-  mean(air[213:244],na.rm=TRUE)
  gMeanWater <- mean(temps,na.rm=TRUE)
  gMedianWater <- median(temps,na.rm=TRUE)
  summerairJA <- mean(air[274:335],na.rm=TRUE)
  winterairJF <- mean(air[93:151],na.rm=TRUE)
  winterairDJ <- mean(air[62:120],na.rm=TRUE)
  gMeanAir <- mean(air,na.rm=TRUE)
  gMedianAir <- median(air,na.rm=TRUE)
  
  metrics <- rbind.data.frame(metrics,
                              cbind.data.frame(year=y,gMeanAir,gMedianAir,gMeanWater,gMedianWater,springairA,summerairJA,winterairDJ,winterairJF), make.row.names= F)
}


#==== Snow Water Equivalent ====
# Load file
sweApr1 <- read.csv("./data/sweApril1.csv")
metricz <- left_join(metrics,sweApr1 %>% dplyr::select(c(year,mean)),by="year") %>% rename(swe=mean) 

