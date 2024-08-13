## this is a "lean" version of CarBayesOM.r, just to prepare data 

#  see also "lean start" line 273 
 
 
 setwd("C:/localH/research/Awal/codeFeb15")


library(spdep)
library(CARBayesdata)
library(sf)
require(maps)
require(sp)
require(maptools)
library(leaflet)
library(dplyr)
require(Matrix)
require(tibble)



library(raster)                                         # this one may be time-consuming
st.sp <- getData('GADM', country='USA', level=1)
st.sp <- subset(st.sp, NAME_1 != "District of Columbia")
st.nb <- poly2nb(st.sp)
st.nb <- subset(st.nb, subset=card(st.nb) > 0)
st.mat <- nb2mat(st.nb, style="B",zero.policy = TRUE)
st.mat[(st.mat>0)] = 1



 ############## manage county level data

#### reading vaccination data

d_county <- read.csv ("data_county_timeseries.csv", header = T)
d_county <- subset(d_county, STATE_NAME != "DC")
d_county <- subset(d_county, STATE_NAME != "GU")
d_county <- subset(d_county, STATE_NAME != "PR")
d_county <- subset(d_county, STATE_NAME != "VI")
d_county <- subset (d_county, CASE_TYPE != "Booster")
d_county <- subset (d_county, CASE_TYPE != "Booster Coverage")
d_county <- subset (d_county, CASE_TYPE != "Complete")
d_county <- subset (d_county, CASE_TYPE != "Partial")
d_county <- subset (d_county, CASE_TYPE != "Partial Coverage")
d_vac_county <- d_county[d_county$DATE > "2021-06-29" &  d_county$DATE < "2021-07-07",] ## complete coverage till 07-07-2021
loving <- data.frame ("TX", 48, "Loving", 48301, "County", "3/27/22", "Complete Coverage", 20.1183431952662, 169,	12,	2022)
names(loving) <- c("STATE_NAME",	"STATE",	"COUNTY_NAME",	"COUNTY",	"GEOFLAG",	"DATE",	"CASE_TYPE",	"CASES",	"POPN",	"WEEK",	"YEAR")
d_vac_county <- rbind(d_vac_county, loving)
d_vac_county <- d_vac_county[order(d_vac_county$STATE_NAME, d_vac_county$COUNTY_NAME),]


######## reading case data

d_case_county <- read.csv ("county_cases.csv", header = T)
d_case_county <- subset(d_case_county, Admin2 != "Unassigned")
d_case_county <- subset(d_case_county, Province_State != "Puerto Rico")
d_case_county <- subset(d_case_county, FIPS != "NA")
d_case_county <- subset(d_case_county,  Province_State != "American Samoa")
d_case_county <- d_case_county[!grepl("Out of", d_case_county$Admin2),]
d_case_county <- subset(d_case_county,  Province_State != "Northern Mariana Islands")
d_case_county <- subset(d_case_county, Admin2 != "Bristol Bay plus Lake and Peninsula")
d_case_county <- subset(d_case_county, Admin2 != "Chugach")
d_case_county <- d_case_county[order(d_case_county$State, d_case_county$Admin2),]

######### merging case and vaccination data
dim(d_case_county)
dim(d_vac_county)
d_case_county_100k <- d_case_county[,-c(1:12)]/d_vac_county$POPN ## cases per 100,000 k
d_case_county_100k$tot_cases <-  d_case_county_100k[,710] - d_case_county_100k[,650] ## cases from 07-01-2021 to 12-31-2021
d_case_vac_county <- data.frame (d_case_county$State, d_case_county$Admin2, 
                                 d_case_county_100k$tot_cases, d_vac_county$CASES,
                                 d_case_county$Lat, d_case_county$Long_)
names(d_case_vac_county)[1] <- "state"
names(d_case_vac_county)[2] <- "county"
names(d_case_vac_county)[3] <- "cases"
names(d_case_vac_county)[4] <- "vaccination"
names(d_case_vac_county)[5] <- "Latitute"
names(d_case_vac_county)[6] <- "Longitude"
head(d_case_vac_county)

######### reading demographic data

d_dem_county <- read.csv ("acs2017_county_data.csv", header = T)
d_dem_county <- subset (d_dem_county, State != "Puerto Rico")
d_dem_county <- subset (d_dem_county, State != "District of Columbia")
d_dem_county <- d_dem_county[order(d_dem_county$state, d_dem_county$County),]
head(d_dem_county)

######## merging case, vaccination and demographic data
dim (d_case_vac_county)
dim(d_dem_county)
covidnew <- cbind(d_case_vac_county, d_dem_county)
head(covidnew)
covidnew <- covidnew[,c(1:6,8,11:44)]
head(covidnew)

######################## 

covidnew1 <- read.csv ("covid_us_county.csv", header = T)
dim(covidnew1)
covidnew1 <- covidnew1[covidnew1$date == "2022-10-19",]
covidnew1 <- subset(covidnew1, state != "Virgin Islands")
covidnew1 <- subset(covidnew1, state != "Puerto Rico")
covidnew1 <- subset(covidnew1, state != "Northern Mariana Islands")
covidnew1 <- subset(covidnew1, state != "Guam")
covidnew1 <- subset(covidnew1, state != "Grand Princess")
covidnew1 <- subset(covidnew1, state != "District of Columbia")
covidnew1 <- subset(covidnew1, state != "Diamond Princess")
covidnew1 <- subset(covidnew1, state != "American Samoa")
covidnew1 <- subset(covidnew1, county != "Unassigned")
covidnew1 <- covidnew1[!grepl("Out of", covidnew1$county),]
covidnew1 <- subset(covidnew1 , fips != "NA")
covidnew1 <- subset(covidnew1, county != "Bristol Bay")
covidnew1 <- subset(covidnew1, county != "Chugach")
covidnew1 <- covidnew1[order(covidnew1$state_code, covidnew1$county),]

new_row <- c(fips=49001, county="Beaver",  state="Utah", lat = 38.35657,
             long = -113.2342, date="2020-03-29", cases = 1, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49005, county="Cache",  state="Utah", lat = 41.72330587,
             long = -111.7443667, date="2020-04-16", cases = 37, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49007, county="Carbon",  state="Utah", lat = 39.64834818,
             long = -110.5872512, date="2020-04-16", cases = 2, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49013, county="Duchesne",  state="Utah", lat = 40.29772254,
             long = -110.425237, date="2020-04-16", cases = 2, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49015, county="Emery",  state="Utah", lat = 38.99617072,
             long = -110.7013958, date="2020-04-16", cases = 3, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49017, county="Garfield",  state="Utah", lat = 37.85447192,
             long = -111.4418764, date="2020-04-16", cases = 1, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49019, county="Grand",  state="Utah", lat = 38.98103848,
             long = -109.5704487, date="2020-04-16", cases = 1, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49021, county="Iron",  state="Utah", lat = 37.8590362,
             long = -113.2897381, date="2020-04-16", cases = 15, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49025, county="Kane",  state="Utah", lat = 37.28507306,
             long = -111.8861752, date="2020-04-16", cases = 3, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49047, county="Uintah",  state="Utah", lat = 40.12491499,
             long = -109.5174415, date="2020-04-16", cases = 6, state_code = "UT", deaths = 0)
covidnew1 <- rbind(covidnew1, new_row)
new_row <- c(fips=49053, county="Washington",  state="Utah", lat = 37.28003525,
             long = -113.504698, date="2020-04-18", cases = 77, state_code = "UT", deaths = 1)
covidnew1 <- rbind(covidnew1, new_row)
dim(covidnew1)


############ deletion of 0 cases

cou <- c("Beaver", "Cache", "Carbon", "Duchesne", "Emery",
         "Garfield", "Grand", "Iron", "Kane", "Uintah", "Washington")
num <- NULL
for (i in 1:11){
  num [i] <- which(covidnew1$county == cou[i] & covidnew1$cases == 0)
}
covidnew1 <- covidnew1[-num, ]
head(covidnew1)
dim(covidnew1)
###################################

covidnew1 <- data.frame(as.numeric(covidnew1$fips), covidnew1$county, covidnew1$state, as.numeric(covidnew1$lat),
                        as.numeric(covidnew1$long), as.Date(covidnew1$date), as.numeric(covidnew1$cases),
                        covidnew1$state_code, as.numeric(covidnew1$deaths))

names(covidnew1)[1] <- "fips"; names(covidnew1)[2] <- "county"; names(covidnew1)[3] <- "state"
names(covidnew1)[4] <- "lat" ; names(covidnew1)[5] <- "long"; names(covidnew1)[6] <- "date"
names(covidnew1)[7] <- "cases" ; names(covidnew1)[8] <- "state_code" ; names(covidnew1)[9] <- "deaths"
str(covidnew1)

covidnew1 <- covidnew1[order(covidnew1$state, covidnew1$county),]

covidnew1 <- data.frame(covidnew, covidnew1$cases, covidnew1$cases/covidnew$TotalPop, 
                        covidnew1$deaths, covidnew1$deaths/covidnew$TotalPop)
names(covidnew1)[42] <- "Cases"
names(covidnew1)[43] <- "cases_per100k"
names(covidnew1)[44] <- "Deaths"
names(covidnew1)[45] <- "deaths_per100k"

# mydata extra
# Alaska = 3, "Copper River" "Hoonah-Angoon" "Yakutat"
# Maryland=1, Baltimore City 
# missouri=1, "St. Louis City"

covidnew1[covidnew1$county=="Copper River",]
covidnew1 <- subset(covidnew1, county != "Copper River")

covidnew1[covidnew1$county=="Hoonah-Angoon",]
covidnew1 <- subset(covidnew1, county != "Hoonah-Angoon")

covidnew1[covidnew1$county=="Yakutat",]
covidnew1 <- subset(covidnew1, county != "Yakutat")

covidnew1[covidnew1$county=="Baltimore City",]
covidnew1 <- subset(covidnew1, county != "Baltimore City")

covidnew1[covidnew1$county=="St. Louis City",]
covidnew1 <- subset(covidnew1, county != "St. Louis City")

write.csv(covidnew1, "new county level final data.csv")
covidnew1 <- read.csv ("new county level final data.csv", header = T)


########## county neighbor

# neighbor data extra
# district of columbia = 1
# Illinois = 1, Indiana=1, michigan = 4, minnesota = 1, newyork=1
# Ohio =1, Wisconsin =2, 


covidnew <- covidnew[order(covidnew$State, covidnew$county),]   # do I need this?

library(raster) # loads shapefile
co.sp <- getData('GADM', country='USA', level=2)  
  #  getData will be removed in a future version of raster. Please use the geodata package instead

# covidnew[covidnew$State=="Missouri",]$county
# co.sp[co.sp$NAME_1=="Missouri",]$NAME_2

co.sp <- subset(co.sp, NAME_1!="District of Columbia")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.14.51_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.15.46_1")

co.sp[co.sp$NAME_2=="Lake Hurron",]
co.sp <- subset(co.sp, GID_2 != "USA.23.44_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.23.45_1")

co.sp[co.sp$NAME_2=="Lake St. Clair",]
co.sp <- subset(co.sp, GID_2 != "USA.23.46_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.23.47_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.24.40_1")

co.sp[co.sp$NAME_2=="Lake Ontario",]
co.sp <- subset(co.sp, GID_2 != "USA.33.25_1")

co.sp[co.sp$NAME_2=="Lake Erie",]
co.sp <- subset(co.sp, GID_2 != "USA.36.44_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.50.34_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.50.35_1")

co.nb <- poly2nb(co.sp)
co.nb <- subset(co.nb, subset=card(co.nb) > 0)   # co.sp was 3136 counties, after this only 3127! 
co.mat <- nb2mat(co.nb, style="B")

  # saving the neighbor matrix 
  
  Cmat = as(co.mat, "CsparseMatrix")
  writeMM(Cmat,file='Cmat.txt')
  
   # Cmat1 = readMM(file='Cmat.txt')   # to read 
   
   
## ================== "lean" START 

library(spdep)
library(CARBayes)
# library(CARBayesdata)
library(sf)
require(maps)
require(sp)
require(maptools)
library(leaflet)
library(dplyr)
require(Matrix)

 setwd("C:/localH/research/Awal/codeFeb15")

covidnew2 <- read.csv ("covidnew1-st1.csv", header = T)   # formerly "new county level final data.csv"
#covidnew1 <- covidnew1[order(covidnew1$State, covidnew1$county),]

 Cmat2 = readMM(file='Cmat-st1.txt') 
 
 # the datasets don't match, so here's some manupulation 
 
 co.sp1 <- subset(co.sp, NAME_1 == "Indiana")
 dim(co.sp1)
 set1 = covidnew2[ covidnew2$state == "IN",]
 c1 = set1$county 
  
  c2 = co.sp1$NAME_2   # list of Counties 

   levels(as.factor(covidnew2$state))
   levels(as.factor(co.sp$NAME_1))
  
  setdiff(c1,c2)  
  setdiff(c2,c1)
  
    # AK: huge mess!   AL:  "Shelby" and "St. Clair" switched 
	# AR:  "St. Francis" out of order
	# AZ: ok  
	# CA: ok   CO: ok   CT: ok   DE: ok
	# FL:  "Santa Rosa"  "Sarasota"     "Seminole"     "St. Johns"    "St. Lucie"  scrambled 
	# GA: ok  HI: ok  IA: ok  ID: ok   
	# IL:  counties between "De Kalb"     "Dupage"      "La Salle"    "Saint Clair" scrambled 
	# IN: "De Kalb"      "Saint Joseph" 
	# KS: 

  ## also: the order of some states are ALSO not the same in the two datasets!  (e.g. try Virginia and Vermont) 
  #   (Alaska & Alabama switched too) 

 
# [1] "AK" "AL" "AR" "AZ" "CA" "CO" "CT" "DE" "FL" "GA" "HI" "IA" "ID" "IL"
#[15] "IN" "KS" "KY" "LA" "MA" "MD" "ME" "MI" "MN" "MO" "MS" "MT" "NC" "ND"
#[29] "NE" "NH" "NJ" "NM" "NV" "NY" "OH" "OK" "OR" "PA" "RI" "SC" "SD" "TN"
#[43] "TX" "UT" "VA" "VT" "WA" "WI" "WV" "WY"
 
 
 # a routine to gather all "bad counties"
 
 st.list =  levels(as.factor(co.sp$NAME_1)) 
 tot.bad = 0
 for (i in 3:length(st.list)){
     st.name = st.list[i]
	 co.sp1 <- subset(co.sp, NAME_1 == st.name)
	 dim(co.sp1)
	 set1 = covidnew2[ covidnew2$State == st.name,]
	 c1 = set1$county 
	 c2 = co.sp1$NAME_2   # list of Counties 
	  print(st.name)
	  print(setdiff(c2,c1))
	  print(setdiff(c1,c2) ) 
	  tot.bad = tot.bad + length(setdiff(c2,c1))
 } 
 
 ## begin synching the two lists 
 
  # purge Alaska 
   co.sp2 = subset(co.sp, NAME_1 != "Alaska")
   covidnew2a = covidnew2[ covidnew2$State != "Alaska",]
   Nrec =  dim(co.sp2)[1]    # 3110 counties
   sublist = read.table("sub-list1.txt", sep=",", header = F, strip.white = T)  
        # 2nd column: county name in "co.sp", 3rd column: county name in "covidnew1"
   
  # main loop
     n.fatal = 0
	 n.fatal2 = 0
     covidnew3 = covidnew2a
  
    for (k in 1:Nrec){
	  rec = covidnew2a[k,]
	  k1 = which((co.sp2$NAME_1 == rec$State) & (co.sp2$NAME_2 == rec$county))
	  if (length(k1) == 0){     # no match 
	     k2 = which((sublist$V1 == rec$State) & (sublist$V3 == rec$county))
		 if (length(k2) == 0){     # no match
		   	n.fatal = n.fatal + 1 
		 }  else {  # if found on the sub lsit
		     new.name = sublist$V2[k2]
	         k3 = which((rec$State == co.sp2$NAME_1) & (new.name == co.sp2$NAME_2))
			 if (length(k2) == 0){     # no match
		   	    n.fatal2 = n.fatal2 + 1 
		      }  else {  #write into new, transforming the name
			     covidnew3[k3,] = covidnew2a[k,]
				 covidnew3$county[k3] = new.name
			  }
	     }
	  }  else {
	    covidnew3[k1,] = covidnew2a[k,]   #write into new 
	  }
	}    
   c(n.fatal, n.fatal2)   

   # for testing
   
   subset(co.sp2, NAME_1 == "Virginia")$NAME_2
   covidnew3[ covidnew3$State == "Virginia",]$county      # these should match 
   
   
    # which((co.sp$NAME_1 == "Alabama") & (co.sp$NAME_2 == "De Kalb"))
   ## Virginia has 2 unmatched counties, need to work manually on that! 
   #    1) remove from co.sp:  "Clifton Forge City",  "Bedford City" 
   #     2) add "Richmond" (city, not county) 3) not sure what to do with "Franklin City", maybe just delete from covidnew2
  
   # write.csv(covidnew3, "covidnew3.csv") 
    
 co.nb2 <- poly2nb(co.sp2)
   list0 = which( card(co.nb2) == 0)
 co.nb2 <- subset(co.nb2, subset=card(co.nb2) > 0) 
 co.mat2 <- nb2mat(co.nb2, style="B")   # here, we cannot allow counties with 0 neighbors

 covidnew3$logc = log(covidnew3$cases_per100k + 0.001)
 covidnew3$vaxx = covidnew3$vaccination/100
 covidnew3$vaxx0 = covidnew3$vaxx - mean(covidnew3$vaxx)
  covidnew3$logd = log(covidnew3$deaths_per100k + 0.001)
 
  covidnew3a = covidnew3[-list0,]    # renamed to 3a (inteferes with East/West stuff)
  
    # now ready to run myownCAR1.r 
  
 
 ## ===========================================================
 
 covidnew1$logc = log(covidnew1$cases_per100k + 0.001)
 covidnew1$vaxx = covidnew1$vaccination/100
 covidnew1$vaxx0 = covidnew1$vaxx - mean(covidnew1$vaxx)
 
(chain_case <- S.CARleroux(logc ~ vaxx0 + state, data=covidnew1,family="gaussian",W= 0+ as.matrix(Cmat),
             prior.tau2 = c(10,1), prior.nu2 = c(10,0.01), rho = 0.95, burnin=1000, n.sample=50000, thin = 10))

 samp = chain_case$samples

 ## =========== for the East/West stuff 

 subsWest = c("WA","OR","CA","AZ","NM","TX", "OK", "NE", "SD", "ND", "MT", "ID", "NV", "UT", "WY", "CO", "KS") 
 ww = which(is.element(covidnew3$state, subsWest))
 cnew.west = covidnew3[ww,]
 
 west <- c("Washington", "Oregon", "California", "Arizona", "New Mexico", "Texas",
          "Oklahoma","Nebraska", "South Dakota", "North Dakota","Montana", "Idaho",
          "Nevada", "Utah", "Wyoming", "Colorado", "Kansas")

  co.sp_west <- subset(co.sp2, NAME_1 %in% west)
  co.nb_west <- poly2nb(co.sp_west)
   list0w = which(card(co.nb_west) == 0)
  co.nb_west <- subset(co.nb_west, subset=card(co.nb_west) > 0)
  co.mat_west <- nb2mat(co.nb_west, style="B")
    cnew.west = cnew.west[-list0w,]   # skip this, if it appears I already skipped that county when updating covidnew3
  dim(co.mat_west)
  dim(cnew.west)
 
   # now ready to run myownCAR1w.r 
 
  cnew.east = covidnew3[-ww,]

 co.sp_east <- subset (co.sp2, !(NAME_1 %in% west) & !(NAME_1 == "Hawaii"))
   cnew.east = cnew.east[cnew.east$state != "HI",]
   co.sp_east$yhat1 = cnew.east$logc
   spplot(co.sp_east,"yhat1")
	
 co.nb_east <- poly2nb(co.sp_east)
   list0e = which(card(co.nb_east) == 0)
 co.nb_east <- subset(co.nb_east, subset=card(co.nb_east) > 0)
 co.mat_east <- nb2mat(co.nb_east, style="B")

 cnew.east = cnew.east[-list0e,]
  dim(co.mat_east)
  dim(cnew.east)

     # now ready to run myownCAR1e.r 
 
  
 
(chain_case <- S.CARleroux(logc ~ vaxx0 + state, data=covidnew1_west,family="gaussian",W= co.mat_east,
             prior.tau2 = c(10,1), prior.nu2 = c(10,0.01), rho = 0.95, burnin=1000, n.sample=50000, thin = 10))


## a "small" dataset to debug my own MCMC routine 

co.spTX <- subset(co.sp, NAME_1 == "Texas")
co.nbTX <- poly2nb(co.spTX)
co.nbTX <- subset(co.nbTX, subset=card(co.nbTX) > 0)
co.matTX <- nb2mat(co.nbTX, style="B")

covidnew1TX <- covidnew1[covidnew1$state == "TX",]

(chain_case <- S.CARleroux(logc ~ vaxx0, data=covidnew1TX, family="gaussian", W= co.matTX,
             prior.tau2 = c(10,1), prior.nu2 = c(10,0.01), rho = 0.95, burnin=1000, n.sample=50000, thin = 10))

# ==============
# quick start for the SP data 

library(spdep)
library(CARBayes)
# library(CARBayesdata)
library(sf)
require(maps)
require(sp)
require(maptools)
library(leaflet)
library(dplyr)
require(Matrix)

library(raster) # loads shapefile
co.sp <- getData('GADM', country='USA', level=2)  
  #  getData will be removed in a future version of raster. Please use the geodata package instead

# covidnew[covidnew$State=="Missouri",]$county
# co.sp[co.sp$NAME_1=="Missouri",]$NAME_2

co.sp <- subset(co.sp, NAME_1!="District of Columbia")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.14.51_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.15.46_1")

co.sp[co.sp$NAME_2=="Lake Hurron",]
co.sp <- subset(co.sp, GID_2 != "USA.23.44_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.23.45_1")

co.sp[co.sp$NAME_2=="Lake St. Clair",]
co.sp <- subset(co.sp, GID_2 != "USA.23.46_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.23.47_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.24.40_1")

co.sp[co.sp$NAME_2=="Lake Ontario",]
co.sp <- subset(co.sp, GID_2 != "USA.33.25_1")

co.sp[co.sp$NAME_2=="Lake Erie",]
co.sp <- subset(co.sp, GID_2 != "USA.36.44_1")

co.sp[co.sp$NAME_2=="Lake Michigan",]
co.sp <- subset(co.sp, GID_2 != "USA.50.34_1")

co.sp[co.sp$NAME_2=="Lake Superior",]
co.sp <- subset(co.sp, GID_2 != "USA.50.35_1")

  co.sp2 = subset(co.sp, NAME_1 != "Alaska")
  covidnew3 = read.csv("covidnew3.csv") 
 
  co.nb2 <- poly2nb(co.sp2)
   list0 = which( card(co.nb2) == 0)
 co.nb2 <- subset(co.nb2, subset=card(co.nb2) > 0) 
 co.mat2 <- nb2mat(co.nb2, style="B")   # here, we cannot allow counties with 0 neighbors

 covidnew3$logc = log(covidnew3$cases_per100k + 0.001)
 covidnew3$vaxx = covidnew3$vaccination/100
 covidnew3$vaxx0 = covidnew3$vaxx - mean(covidnew3$vaxx)
  covidnew3$logd = log(covidnew3$deaths_per100k + 0.001)
 
  covidnew3a = covidnew3[-list0,]    # renamed to 3a (inteferes with East/West stuff)
  
   dim(co.mat2) 
   dim(covidnew3a)
   
    # now ready for myownCAR1.r, or additionally do East/West stuff
	
	## EZ plotting 
	
	co.sp_west$yhat1 = cnew.west$vaxx
    spplot(co.sp_west,"yhat1")
	
	
   
  ## ============= NEW data prep: more months for cases and vacc., other covariates 
  
  #### reading vaccination data

d_county <- read.csv ("data_county_timeseries.csv", header = T)
d_county <- subset(d_county, STATE_NAME != "DC")
d_county <- subset(d_county, STATE_NAME != "GU")
d_county <- subset(d_county, STATE_NAME != "PR")
d_county <- subset(d_county, STATE_NAME != "VI")
d_county <- subset (d_county, CASE_TYPE != "Booster")
d_county <- subset (d_county, CASE_TYPE != "Booster Coverage")
d_county <- subset (d_county, CASE_TYPE != "Complete")
d_county <- subset (d_county, CASE_TYPE != "Partial")
d_county <- subset (d_county, CASE_TYPE != "Partial Coverage")
vac.Jul <- d_county[d_county$DATE > "2021-06-29" &  d_county$DATE < "2021-07-07",] ## complete coverage till 07-07-2021
 dd = which(duplicated(vac.Jul$COUNTY))

vac.Aug <- d_county[d_county$DATE > "2021-07-29" &  d_county$DATE < "2021-08-07",]
 dd = which(duplicated(vac.Aug$COUNTY))
 vac.Aug[2932,]
 vac.Aug = vac.Aug[!((vac.Aug$COUNTY == 51678) & ( vac.Aug$WEEK == 31)),]
 
vac.Sep <- d_county[d_county$DATE > "2021-08-29" &  d_county$DATE < "2021-09-07",]
 dd = which(duplicated(vac.Sep$COUNTY))
 
  vac.Oct <- d_county[d_county$DATE > "2021-09-29" &  d_county$DATE < "2021-10-07",]
 dd = which(duplicated(vac.Oct$COUNTY))
 
  # vac.Nov <- d_county[d_county$DATE > "2021-10-29" &  d_county$DATE < "2021-11-07",]
 # dd = which(duplicated(vac.Nov$COUNTY))
   # badlist = vac.Nov$COUNTY[duplicated(vac.Nov$COUNTY)]
   # vac.Nov[vac.Nov$COUNTY %in% badlist,]
   # vac.Nov = vac.Nov[!((vac.Nov$COUNTY %in% badlist) & ( vac.Nov$WEEK == 44)),]
   # # one missing value 
    # setdiff(vac.Oct$COUNTY, vac.Nov$COUNTY)
	
  	
	 vac.Nov <- d_county[d_county$WEEK == 44,]  # this is way easier! 
    dd = which(duplicated(vac.Nov$COUNTY))
	

	vac.Dec <- d_county[d_county$WEEK == 48,]  # this is way easier! 
    dd = which(duplicated(vac.Dec$COUNTY))
	
	JulCASES = vac.Jul$CASES
	vac.Jul$CASES = NULL	
	vac.All = cbind(vac.Jul, JulCASES, vac.Aug$CASES, vac.Sep$CASES, vac.Oct$CASES, vac.Nov$CASES, vac.Dec$CASES) 
	names(vac.All)[11:16] = c("vacJul", "vacAug","vacSep","vacOct","vacNov","vacDec")
	
	nli = 40
	rndCty = sample(1:3140,nli,replace = F) 
	
	plot(1:6,vac.All[rndCty[1], 11:16],ylim=c(0,100))
	for (i in 2:nli){
	   lines(1:6,vac.All[rndCty[i], 11:16],col= rgb(i/nli,0.5,1 - i/nli))
	}
	
	
	## for the case some counties have a DECREASING vacc %? 
	
	# E.g.   >  vac.All[2529,]
       # STATE_NAME STATE COUNTY_NAME COUNTY      GEOFLAG       DATE
       #         TX    48    Atascosa  48013 Interpolated 2021-07-04
       #                CASE_TYPE  POPN WEEK YEAR   vacJul   vacAug   vacSep   vacOct
       # 919633 Complete Coverage 51153   26 2021 54.13816 51.66162 48.56595 46.08941
       #  vacNov   vacDec
       #  44.93774 46.14588
  
 
 loving <- data.frame("TX", 48, "Loving", 48301, "County", "3/27/22", "Complete Coverage", 169,	12,	2022, 20.118, 20.118, 20.118, 20.118, 20.118, 20.118)
 names(loving) = names(vac.All)
vac.All  <- rbind(vac.All, loving)
vac.All <- vac.All[order(vac.All$STATE_NAME, vac.All$COUNTY_NAME),]
 
 

######## reading case data: these two are copied one for one

d_case_county <- read.csv ("county_cases.csv", header = T)
d_case_county <- subset(d_case_county, Admin2 != "Unassigned")
d_case_county <- subset(d_case_county, Province_State != "Puerto Rico")
d_case_county <- subset(d_case_county, FIPS != "NA")
d_case_county <- subset(d_case_county,  Province_State != "American Samoa")
d_case_county <- d_case_county[!grepl("Out of", d_case_county$Admin2),]
d_case_county <- subset(d_case_county,  Province_State != "Northern Mariana Islands")
d_case_county <- subset(d_case_county, Admin2 != "Bristol Bay plus Lake and Peninsula")
d_case_county <- subset(d_case_county, Admin2 != "Chugach")
d_case_county <- d_case_county[order(d_case_county$State, d_case_county$Admin2),]
  
   ### 
   dim(d_case_county)
   dim(d_vac_county)
d_case_county_100k <- d_case_county[,-c(1:12)]/d_vac_county$POPN ## cases per 100,000 k

  namc = names(d_case_county)
  # namc[691]   # Nov.30     #  namc[722]    # Dec.31
  
  case.pAug = d_case_county_100k[,600] - d_case_county_100k[,569]
  case.pSep = d_case_county_100k[,630] - d_case_county_100k[,600]
  case.pOct = d_case_county_100k[,661] - d_case_county_100k[,630]
  case.pNov = d_case_county_100k[,691] - d_case_county_100k[,661]
  case.pDec = d_case_county_100k[,722] - d_case_county_100k[,661]
  case.pJan = d_case_county_100k[,753] - d_case_county_100k[,722]  
  
  #d_case_county_100k$tot_cases <-  d_case_county_100k[,710] - d_case_county_100k[,650] ## cases from 07-01-2021 to 12-31-2021 ??-OM
d_case_vac_county <- data.frame (d_case_county$State, d_case_county$Admin2, d_case_county$Lat, d_case_county$Long_,
                                 case.pAug, case.pSep, case.pOct, case.pNov, case.pDec, case.pJan)
								 								 						 
names(d_case_vac_county)[1:4] <- c("state", "county","Latitude", "Longitude")

d_case_vac_county = cbind(d_case_vac_county, vac.All[,11:16])

head(d_case_vac_county)

######### reading demographic data

d_dem_county <- read.csv ("acs2017_county_data.csv", header = T)
d_dem_county <- subset (d_dem_county, State != "Puerto Rico")
d_dem_county <- subset (d_dem_county, State != "District of Columbia")
d_dem_county <- d_dem_county[order(d_dem_county$state, d_dem_county$County),]
head(d_dem_county)

######## merging case, vaccination and demographic data
dim (d_case_vac_county)
dim(d_dem_county)
covidnew <- cbind(d_case_vac_county, d_dem_county)
head(covidnew)
covidnew <- covidnew[,-(19:20)]
head(covidnew)

###### add death data 

if (0==1){

covidnew1 <- read.csv ("covid_us_county.csv", header = T)   # this one contains cumulative cases & deaths
dim(covidnew1)
dethAug = covidnew1$deaths[covidnew1$date == "2021-08-31"] - covidnew1$deaths[covidnew1$date == "2021-07-31"]
dethSep = covidnew1$deaths[covidnew1$date == "2021-09-30"] - covidnew1$deaths[covidnew1$date == "2021-08-31"]
dethOct = covidnew1$deaths[covidnew1$date == "2021-10-31"] - covidnew1$deaths[covidnew1$date == "2021-09-30"]
dethNov = covidnew1$deaths[covidnew1$date == "2021-11-30"] - covidnew1$deaths[covidnew1$date == "2021-10-31"]
dethDec = covidnew1$deaths[covidnew1$date == "2021-12-31"] - covidnew1$deaths[covidnew1$date == "2021-11-30"]
dethJan = covidnew1$deaths[covidnew1$date == "2022-01-31"] - covidnew1$deaths[covidnew1$date == "2021-12-31"]

covidnew1 = covidnew1[covidnew1$date == "2021-08-31",c(1:5,8)]
covidnew1a = cbind(covidnew1, dethAug, dethSep, dethOct, dethNov, dethDec, dethJan)

  drop.list = c("Virgin Islands", "Puerto Rico", "Northern Mariana Islands", "Guam", "Grand Princess", "District of Columbia",
     "Diamond Princess"," American Samoa") 
   covidnew1 <- subset(covidnew1, !is.element(state, drop.list))  
   covidnew1 <- subset(covidnew1, county != "Unassigned")
   covidnew1 <- covidnew1[!grepl("Out of", covidnew1$county),]
    covidnew1 <- subset(covidnew1 , fips != "NA")
    covidnew1 <- subset(covidnew1, county != "Bristol Bay")
    covidnew1 <- subset(covidnew1, county != "Chugach")
  covidnew1 <- covidnew1[order(covidnew1$state_code, covidnew1$county),]
  
  # cleanup of negative county resutls 
  
  covidnew1 <- read.csv ("covid_us_county.csv", header = T) 
  c1 = covidnew1[covidnew1$fips == 6085,]    #  Santa Clara CA, one of the worst! 
  
   ## Dukes and Nantucket counties, MA, are mashed together, for some reason:   starts at row   1271539 
   #    starts at row   1308613  etc... 
   
   # simplest option: delete the NA's 
   
   covidnew1 = covidnew1[!is.na(covidnew1$fips),]
   cnew1 = covidnew1
    c1 = cnew1[cnew1$fips == 5007,] 
   
   c1d = as.Date(c1$date, format = '%Y-%m-%d')
   c1di = as.integer(c1d)
   c1 = c1[order(c1di),]    # sorted by date 
   plot(c1$deaths) 
      n = dim(c1)[1]
	  c1inc = c1$deaths[2:n] - c1$deaths[1:(n-1)]
	}

	
   ## trying the new dataset now! 
   
   cnew2 = read.csv("time_series_covid19_deaths_US.csv") 
     c1 = cnew2[216,13:1000]
     plot(as.numeric(c1))       #the same problem persists! 

	cnew2 = cnew2[!is.na(cnew2$FIPS),]
	cnew2 = cnew2[cnew2$Admin2 != "Unassigned",] 
 	cnew2 = cnew2[!grepl("Out of", cnew2$Admin2),]
	cnew2 = cnew2[cnew2$Province_State != "Puerto Rico",]
			 
	# monitoring large negatives 

	 fips = cnew2$FIPS
	nfips = length(fips)
	endT = 1155
	T = endT - 12
	negsum = numeric(nfips)
	bigj = numeric(nfips);   maxd = numeric(nfips)
	for (k in 1:nfips){
	   fk = fips[k]
	   c1 = cnew2[fips == fk,] 
	   deth = as.numeric(c1[13:(T+12)])
	   dinc = deth[2:T] - deth[1:(T-1)]
	   negsum[k] = sum(dinc[dinc < 0])
	   bigj[k] = max(dinc)
	   maxd[k] = max(deth)
	}
	
	 n100 = which(negsum < -100)            # Large negs list: 235, 237, many FL counties...  ==> bad FIPS are 6081, 6085 
	  for (i in 1:12){
	    k =  n100[i+12]
	    c1 = as.numeric(cnew2[k,13:endT])
	    plot(c1, main = cnew2$Combined_Key[k])	 
	 }
	                        
     ratio = bigj/(maxd+1)
     r03 = which(ratio > 0.3)    # Large jumps:  154 163  334-... (all of Florida!)  "Out of ..."  904 
	                              #  968  982 1001 1517 1535 1543 1557 1560 1563 1572 1579 1583 (most of Missouri) 
	                              # most of Nebraska, 1832 
	 par(mfrow=c(4,3))            # ==>  bad FIPS are  5073, 5091, 19185, 20117, 20143, 20179, 34041
     par(mar = c(2,2,4,2))
	 for (i in 1:12){
	    k = r03[i+70]
	    c1 = as.numeric(cnew2[k,13:endT])
	    plot(c1, main = cnew2$Combined_Key[k])	 
	 }
	 
	# now massive purge
	cnew2 = cnew2[cnew2$Province_State != "Florida",]
	cnew2 = cnew2[cnew2$Province_State != "Nebraska",]
	cnew2 = cnew2[cnew2$Province_State != "Missouri",]
	cnew2 = cnew2[cnew2$Province_State != "Alaska",]
	cnew2 = cnew2[cnew2$Province_State != "Hawaii",]
	cnew2 = cnew2[cnew2$Admin2 != "",] 
	cnew2 = cnew2[cnew2$Admin2 != "District of Columbia",]	
	
	cnew1 = covidnew
	cnew1 = cnew1[cnew1$State != "Florida",]
	cnew1 = cnew1[cnew1$State != "Nebraska",]
	cnew1 = cnew1[cnew1$State != "Missouri",]
	cnew1 = cnew1[cnew1$State != "Alaska",]
	cnew1 = cnew1[cnew1$State != "Hawaii",]
	
	sdf1 = setdiff(cnew2$Admin2, cnew1$county)
	sdf2 = setdiff(cnew1$county, cnew2$Admin2)

	 badlist = c(6081, 6085, 5073, 5091, 19185, 20117, 20143, 20179, 34041) 
	  par(mfrow=c(4,3))   
	 for (k in 1:length(badlist)){
	   fk = badlist[k]
	   cnty = which(cnew2$FIPS == fk)
	   c1 = cnew2[cnty,] 
	   deth = as.numeric(c1[13:(T+12)])
	   dinc = deth[2:T] - deth[1:(T-1)]
	   plot(deth,main = cnew2$Combined_Key[cnty])
	}
	
	cnty = which(cnew2$FIPS == 34041)   # fix Warren co. NJ
	cnew2[cnty,771:774] = 333

	# ................ more fixin here
	
  # rescale Santa Clara 
  
      cnty = which(cnew2$FIPS == 6085)
	  dseq = as.numeric(cnew2[cnty, 15:endT])
	   c1 = as.numeric(cnew2[cnty,550:650])
	   plot(c1)   # obs. up to #578 needs to be rescaled 
	   fac = cnew2[cnty,579]/cnew2[cnty,578]
	   # bcup = as.numeric(cnew2[cnty, 15:endT])
	   cnew2[cnty,15:578] = round(cnew2[cnty,15:578]*fac)
	   plot(as.numeric(cnew2[cnty,15:endT))
	   
	   
	 # a more automatic routine? 
          #  cnty = which(cnew2$FIPS == 5073)    cnty = which(cnew2$FIPS == 5091)
	  c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	  dinc = deth[2:T] - deth[1:(T-1)]
      jump = which.max(abs(dinc)) 
      fac = deth[jump+1]/deth[jump]	  
	  newdeth = deth
	  newdeth[1:jump] = round(deth[1:jump]*fac)
	  plot(newdeth) 
	  
	  cnew2[cnty,13:endT] = newdeth
	  plot(as.numeric(cnew2[cnty,15:endT]))
	  
	  # San Mateo CA 
 
      cnty = which(cnew2$FIPS == 6081)
	  c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	  deth[500:510]
	  cnew2[cnty,513:518] = 500 
	  
	  # ignore all the others, their jumps occur earlier 	 
	  
	  ## now remove all negative increments 
	  
	  nc = dim(cnew2)[1]	 
	
	  par(mfrow=c(5,4))
	  par(mar = c(2,2,4,2))
	  for (k in 1:20){
			k1 = k + 2179
			bcup = as.numeric(cnew2[k1,13:endT])
			bcup = c(numeric(12), bcup)
			bcup0 = bcup 		# tracking back the data until hit the negative increment 
			prev = bcup[endT]  #   cnew2[k,endT]
			for (i in (endT-1):13){
			   if (bcup[i] > prev){ 
				 bcup[i] = prev
			   } else {
				 prev = bcup[i]
			   }
			}	  
			 plot(bcup0, main=k1, col="yellow")
			 lines(bcup, col="blue")
	  }
		  
	  
	     ## more bad counties:   1, 100, 111, 123, 132, 148, 169, 209, 637, 787, 1882, 2149, 2150, 2163-2206 many sus counties, 
	        # 2639, 2705, 2710, 2714, 2719
	   
    badlist2 = c(1, 100, 111, 123, 132, 148, 169, 209, 637, 787, 1882, 2149, 2150, 2163, 2206, 2639, 2705, 2710, 2714, 2719)
	# corresponding FIPS 
	badFIPS2 = cnew2$FIPS[badlist2] 
   # 1001  5035  5057  5081  5099  5133  6023  6103 18083 20001 40121 47061 47063 47089 47175 51830 54045 54055 54063 54073
   
   

	  par(mfrow=c(5,4))   
	 for (k in 1:length(badlist2)){
	   cnty = badlist2[k]
	   c1 = cnew2[cnty,] 
	   deth = as.numeric(c1[13:(T+12)])
	   dinc = deth[2:T] - deth[1:(T-1)]
	   plot(deth,main = cnew2$Combined_Key[cnty])
	}
	
	# some isolated low values    cnty = badlist2[9]

	for (k in 17:19){
      cnty = badlist2[k]	
	  c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	  dinc = deth[2:T] - deth[1:(T-1)]
      jump = which.min(dinc)
	  deth[jump+1] = deth[jump]
	  # plot(deth)
	  cnew2[cnty,13:endT] = deth               # bcup2 = cnew2[cnty,]
	  plot(as.numeric(cnew2[cnty,15:endT]))       #    write.csv(cnew2, "cnew2_bcup.csv", row.names=F)
	 } 
	 
	 cnty = badlist2[20]   # this is a tough one 
	 c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	 plot(deth[1110:T])
	 cnew2[cnty,1133] = 39	 
	  
	 # some standard jumps, a lot of them in Tennessee
	 
	 for (k in 1:8){
	  cnty = badlist2[k]	
	  c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	  dinc = deth[2:T] - deth[1:(T-1)]
      jump = which.max(abs(dinc)) 
      fac = deth[jump+1]/deth[jump]	  
	  newdeth = deth
	  newdeth[1:jump] = round(deth[1:jump]*fac)
	  # plot(newdeth) 	  
	  cnew2[cnty,13:endT] = newdeth
	  plot(as.numeric(cnew2[cnty,15:endT]))
	  }
	  
	  which(cnew2$Province_State == "Tennessee")   # all of jumps on the day 700 (Dec. 22)
      badT = c(2121, 2122, 2124, 2125, 2131, 2133, 2135, 2140, 2144, 2147, 2151, 2155, 2158, 2170, 2171, 2174, 2179, 2180, 2181, 2187, 2188, 2190, 2194, 2195, 2202, 2203)       
	  
	  for (k in 14:26){
	  cnty = badT[k]	
	  c1 = cnew2[cnty,]
	  deth = as.numeric(c1[13:(T+12)])
	  dinc = deth[2:T] - deth[1:(T-1)]
      jump = which.max(abs(dinc)) 
      fac = deth[jump+1]/deth[jump]	  
	  newdeth = deth
	  newdeth[1:jump] = round(deth[1:jump]*fac)
	  # plot(newdeth) 	  
	  cnew2[cnty,13:endT] = newdeth
	  plot(as.numeric(cnew2[cnty,15:endT]))
	  }

	#    write.csv(cnew2, "cnew2_bcup.csv", row.names=F)
	
	 # now do the "back tracking" for all the counties 
	 
	 cnew2 = read.csv("cnew2_bcup.csv")
	 
	 nc = dim(cnew2)[1]	 
	
	  par(mfrow=c(5,4))
	  par(mar = c(2,2,4,2))
	  for (k in 1:nc){
			k1 = k 
			bcup = as.numeric(cnew2[k1,13:endT])
			bcup = c(numeric(12), bcup)
			bcup0 = bcup 		# tracking back the data until hit the negative increment 
			prev = bcup[endT]  #   cnew2[k,endT]
			for (i in (endT-1):13){
			   if (bcup[i] > prev){ 
				 bcup[i] = prev
			   } else {
				 prev = bcup[i]
			   }
			}	  
			 plot(bcup0, main=k1, col="yellow")
			 lines(bcup, col="blue")
			 cnew2[k1,13:endT] = bcup[13:endT]   
	  }
		 
	 	# for (k in 1:20){
	  # c1 = cnew2[k,]
	  # deth = as.numeric(c1[13:(T+12)])
	  # dinc = deth[2:T] - deth[1:(T-1)]
      # minc = min(dinc)
	  # print(minc)
	 # }
	 
	#    write.csv(cnew2, "cnew2_final.csv", row.names=F)
	
	dethAug = cnew2[ ,600] - cnew2[,569]
	dethSep = cnew2[ ,630] - cnew2[,600]
	dethOct = cnew2[ ,661] - cnew2[,630]
    dethNov = cnew2[ ,691] - cnew2[,661]
	dethDec = cnew2[ ,722] - cnew2[,691]
	dethJan = cnew2[ ,753] - cnew2[,722]
	dethFeb = cnew2[ ,781] - cnew2[,753]
	
	dethAll = cbind(dethAug, dethSep, dethOct, dethNov, dethDec, dethJan, dethFeb)
	dethPer = cbind(dethAug/cnew2[,12], dethSep/cnew2[,12], dethOct/cnew2[,12], dethNov/cnew2[,12], dethDec/cnew2[,12], dethJan/cnew2[,12], dethFeb/cnew2[,12])*1e5    # deaths per 100,000
	
	monthList = c("dethAug", "dethSep", "dethOct", "dethNov", "dethDec", "dethJan", "dethFeb")
	colnames(dethPer) = monthList
	pairs(dethPer)
	
	kk = apply(dethAll,2,sum)
	
	plot(kk, type="o", xaxt = "n", ylim=c(0,30))
		axis(1, at= 1:7, labels= monthList)
	
	 plot(cnew1$TotalPop, cnew2$Population, log="xy")    # seeing some mismatch! 
	   text(cnew1$TotalPop, cnew2$Population, 1:nc)
	 
	   cnew1 = cnew1[order(cnew1$State,cnew1$county),]	# that cures most of it, but not all    
	   cnew2 = cnew2[order(cnew2$Province_State,cnew2$Admin2),]	
	
	    cnew1$TotalPop[2554] = 56277;  cnew1$TotalPop[2555] =  8334
		cnew2$Population[1038] = 432493; cnew2$Population[1037] = 31368
	    cnew1$TotalPop[2548] = 1142004; cnew1$TotalPop[2549] = 23580
	    cnew1$TotalPop[2615] = 8873;  cnew1$TotalPop[2616] = 220892
	
	
	 # rge = 2401:2600
	 # plot(cnew1$TotalPop[rge], log="y", pch = 16, col="orange") 
	   # lines(cnew2$Population[rge], col="blue", lwd=2) 
	
        # 
	   dethPer = cbind(dethAug/cnew2[,12], dethSep/cnew2[,12], dethOct/cnew2[,12], dethNov/cnew2[,12], dethDec/cnew2[,12], dethJan/cnew2[,12], dethFeb/cnew2[,12])*1e5     # recalculate
		colnames(dethPer) = monthList
		
	cnew1a = cbind(cnew1[,1:16], dethPer, cnew1[,17:52])	   	 
	
	
	vaccdeth = cbind(cnew1[,11:16], dethPer)
	# pairs(vaccdeth) 
	round(cor(vaccdeth),3) 
	
	totDeth = apply(dethPer,1,sum)
	round(cor(cnew1[,11:16], totDeth),3)
	
	vaccase = cbind(cnew1[,11:16], cnew1[,5:10])    # might also fix CasePer100K since some changes were made to the popul.
	round(cor(vaccase),3) 
	
	  ## forgot! add death numbers to *cnew1*
	
	
	
	## synch with the sp list
	
	insert.rc = function(mat, pos){   # insert a row and column into the matrix at "pos"
	   n = dim(mat)[2]
	   if (pos > n){
	      return(0)
		 } else{ 
		   mat1 = cbind(mat,rep(0,n))
		   mat1 = rbind(mat1,rep(0,n+1))		   
		   mat1[pos,1:(n+1)] = 0
		   mat1[1:(n+1),pos] = 0
		   if (pos > 1){
		     mat1[1:(pos-1),(pos+1):(n+1)] = mat[(1:pos-1),pos:n]		   
		     mat1[(pos+1):(n+1), 1:(pos-1)] = mat[pos:n,(1:pos-1)]
		   }	 
		   mat1[(pos+1):(n+1),(pos+1):(n+1)] = mat[pos:n,pos:n]
	       return(mat1)
	     }
	}
	
	
	co.sp3 = co.sp
	dim(co.sp3)
	
	co.sp3 <- subset(co.sp3, NAME_1!="District of Columbia")
	co.sp3 <- subset(co.sp3, NAME_1!="Florida")
	co.sp3 <- subset(co.sp3, NAME_1!="Alaska")
	co.sp3 <- subset(co.sp3, NAME_1!="Hawaii")
	co.sp3 <- subset(co.sp3, NAME_1!="Missouri")
	co.sp3 <- subset(co.sp3, NAME_1!="Nebraska")
	
	dim(co.sp3)
    sdf1 = setdiff(cnew2$Admin2 ,co.sp3$NAME_2)   # "Baltimore City" is missing, will have to add it to co.sp3 manually
    sdf2 = setdiff(co.sp3$NAME_2,cnew2$Admin2)
	
	sp.co = co.sp3$NAME_2
	sp.st = co.sp3$NAME_1
	
	co.nb3 <- poly2nb(co.sp3)
	 nonhb =  which(card(co.nb3) == 0)
    co.nb3 <- subset(co.nb3, subset=card(co.nb3) > 0)
    co.mat3 <- nb2mat(co.nb3, style="B")
	
	dim(co.mat3)  # 2827 2827        4 others missing, add manually?
    dim(cnew2)    # 2832 1155

	    nonhb   # 1118 1124 1560 2670
	
	 sp.co[nonhb]   # "Dukes"     "Nantucket" "Richmond"  "San Juan" 
     sp.st[nonhb]   # "Massachusetts" "Massachusetts" "New York"   "Washington"   

	   ## need to fix: 1) remove Dukes & Nantucket, San Juan
		# 2) add manually the neighbors to Richmond co.NY (aka Staten Island) 
	   #  	New Jersey    Hudson County — north and northeast     Union County — northwest     
	   #  	Middlesex County — west and southwest   Monmouth County — south,   New York     Kings County — east    New York County — northeast
	   
	   sp.co[1090:1095] # "Washington"   "York"   "Allegany"     "Anne Arundel" "Baltimore"  "Calvert" 
	   # Baltimore City is surrounded by Baltimore County, but is politically independent of it. It is bordered by Anne Arundel County
	   sp.co3 = c(sp.co[1:1094],"Baltimore City", sp.co[1095:2831])
	   sp.st3 = c(sp.st[1:1094],"Maryland", sp.st[1095:2831])
	   co.mat3a = insert.rc(co.mat3, 1095)     
	   which(co.sp3$NAME_2 == "Baltimore");   # which(co.sp3$NAME_2 == "Anne Arundel")
	   which(sp.co3 == "Baltimore")
	   co.mat3a[1094,1095] = 1;  co.mat3a[1095,1094] = 1
	   co.mat3a[1093,1095] = 1;  co.mat3a[1095,1093] = 1 
	   
	   # insert Richmond co, NY
	   
	    which((sp.co3 == "Richmond")&(sp.st3 == "New York"))   # 1561
        co.mat3a = insert.rc(co.mat3a, 1559)
		 which(sp.co3 == "Kings")   #  173 1542     I will only add Kings county 
		 co.mat3a[1559,1540] = 1;  co.mat3a[1540, 1559] = 1
	   
	    sp.co3 = sp.co3[-c(1119, 1125, 2671)]   # remove Dukes & Nantucket, San Juan
		sp.st3 = sp.st3[-c(1119, 1125, 2671)] 
		
		## Next, synch with cnew1, cnew2 
				
		> setdiff(co.sp3$NAME_2, cnew1[,2])
 [1] "De Kalb"                "Saint Clair"            "Saint Francis"          "Dupage"                
 [5] "Saint Joseph"           "Saint Bernard"          "Saint Charles"          "Saint Helena"          
 [9] "Saint James"            "Saint John the Baptist" "Saint Landry"           "Saint Martin"          
[13] "Saint Mary"             "Saint Tammany"          "Saint Mary's"           "Saint Louis"           
[17] "Desoto"                 "Debaca"                 "Saint Lawrence"         "Lamoure"               
[21] "Mc Kean"                "Shannon"                "Dewitt"                 "Bedford City"          
[25] "Clifton Forge City"     "Saint Croix"           
> setdiff( cnew1[,2],co.sp3$NAME_2)
 [1] "St. Clair"            "St. Francis"          "DuPage"               "LaSalle"             
 [5] "St. Joseph"           "St. Bernard"          "St. Charles"          "St. Helena"          
 [9] "St. James"            "St. John the Baptist" "St. Landry"           "St. Martin"          
[13] "St. Mary"             "St. Tammany"          "Baltimore City"       "St. Mary's"          
[17] "St. Louis"            "DeSoto"               "De Baca"              "St. Lawrence"        
[21] "LaMoure"              "McKean"               "Oglala Lakota"        "DeWitt"              
[25] "Franklin City"        "Richmond City"        "St. Croix"           
> 
    #add to co.sp3:  "Baltimore City"    "Franklin City"        "Richmond City"  
	#remove from co.sp3: "De Kalb", "Bedford City",  "Clifton Forge City"; 
	   #   right now ad-hoc replaced two virginia counties by 2 others
	# remove  "Dukes"  and "Nantucket" from cnew1
	# sort out De Kalbs: 
  > which(cnew1$county == "DeKalb")  [1]   25  334  513  612 2139
  > which(sp.co3 == "DeKalb")  [1]  334 2137    > which(sp.co3 == "De Kalb") [1]  25 512 610
	
	
	 # setdiff(sp.co3, cnew1[,2])
		
	## repeat the substitution routine. Do I change the full object sp.co3? Or just the matrix (co.mat3)? 	
		
	 Nrec =  dim(co.sp3)[1]    # 2831 counties
    sublist = read.table("sub-list1.txt", sep=",", header = F, strip.white = TRUE)  
        # 2nd column: county name in "co.sp", 3rd column: county name in "covidnew1"
   
  # main loop
     n.fatal = 0
	 n.fatal2 = 0
     cnew3 = cnew1a
  
    for (k in 1:Nrec){
	  rec = cnew1a[k,]
	  k1 = which((co.sp3$NAME_1 == rec$State) & (co.sp3$NAME_2 == rec$county))
	  if (length(k1) == 0){     # no match 
	     k2 = which((sublist$V1 == rec$State) & (sublist$V3 == rec$county))
		 if (length(k2) == 0){     # no match
		   	n.fatal = n.fatal + 1 
		 }  else {  # if found on the sub lsit
		     new.name = sublist$V2[k2]
	         k3 = which((rec$State == co.sp3$NAME_1) & (new.name == co.sp3$NAME_2))
			 if (length(k2) == 0){     # no match
		   	    n.fatal2 = n.fatal2 + 1 
		      }  else {  #write into new, transforming the name
			     cnew3[k3,] = cnew1a[k,]
				 cnew3$county[k3] = new.name
			  }
	     }
	  }  else {
	    cnew3[k1,] = cnew1a[k,]   #write into new 
	  }
	}    
   c(n.fatal, n.fatal2)   # one is "fatal" it's missing Baltimore City
		
		
	 ## check on Lat/Long 	
	 
	 geofit = numeric(Nrec)
	 for (k in 1:Nrec){
	   ee = extent(co.sp3[k,])
       geofit[k] = (ee[1] < cnew3[k,4]) & (ee[2] > cnew3[k,4]) &  (ee[3] < cnew3[k,3]) & (ee[4] > cnew3[k,3])	 
	 }
	 plot(geofit)
	 ww = which(geofit == FALSE)
	 co.sp3[ww,]
	 co.sp3$NAME_2[2820:Nrec]
	 cnew3[2820:Nrec,1:6]      # for some reason, last county got doubled 
	 cnew3[Nrec,] = cnew3[Nrec+1,]
	 
     cnew3 = cnew3[-(Nrec+1),]
	 write.csv(cnew3, "cnew3_woBaltimore.csv", row.names = FALSE)
	 
	 
		# still to do: remove Dukes & Nantucket; add Baltimore City to the matrix co.mat3; fix those 2 Virginia "city" counties.
			
	    # now, repeat the code from earlier (~ line 1120) to make all these adjustments 
		
	sp.co = co.sp3$NAME_2
	sp.st = co.sp3$NAME_1
	
	co.nb3 <- poly2nb(co.sp3)
	 nonhb =  which(card(co.nb3) == 0)
    co.nb3 <- subset(co.nb3, subset=card(co.nb3) > 0)
    co.mat3 <- nb2mat(co.nb3, style="B") 
		
     sp.co[1090:1095] # "Washington"   "York"   "Allegany"     "Anne Arundel" "Baltimore"  "Calvert" 
	   # Baltimore City is surrounded by Baltimore County, but is politically independent of it. It is bordered by Anne Arundel County
	   sp.co3 = c(sp.co[1:1094],"Baltimore City", sp.co[1095:2831])
	   sp.st3 = c(sp.st[1:1094],"Maryland", sp.st[1095:2831])
	   co.mat3a = insert.rc(co.mat3, 1095)     
	   which(co.sp3$NAME_2 == "Baltimore");   # which(co.sp3$NAME_2 == "Anne Arundel")
	   which(sp.co3 == "Baltimore")
	   co.mat3a[1094,1095] = 1;  co.mat3a[1095,1094] = 1
	   co.mat3a[1093,1095] = 1;  co.mat3a[1095,1093] = 1 
	   
	   # insert Richmond co, NY
	   
	    which((sp.co3 == "Richmond")&(sp.st3 == "New York"))   # 1561
        co.mat3a = insert.rc(co.mat3a, 1559)
		 which(sp.co3 == "Kings")   #  173 1542     I will only add Kings county 
		 co.mat3a[1559,1540] = 1;  co.mat3a[1540, 1559] = 1
	   
	    sp.co3 = sp.co3[-c(1119, 1125, 2671)]   # remove Dukes & Nantucket, San Juan
		sp.st3 = sp.st3[-c(1119, 1125, 2671)] 		
		
		          # also remove those from cnew3
		cnew3 = dplyr::add_row(cnew3,cnew1a[1095,], .before = 1095) 
		cnew3 = cnew3[-c(1119, 1125, 2671),]
		
		   # ok now we are synched through cnew3 -- co.mat3a to run analyses
		   # write.csv(cnew3, "cnew3.csv", row.names = FALSE)
		   Cmat3a = as(co.mat3a, "CsparseMatrix")
           writeMM(Cmat3a,file='Cmat3a.txt')
		   
		 totDeth = cnew3$dethAug + cnew3$dethSep + cnew3$dethOct + cnew3$dethNov + cnew3$dethDec + cnew3$dethJan + cnew3$dethFeb
		 cnew3$logd = log(totDeth + 0.1)
		 vaxx = cnew3$vacSep/100
		 cnew3$vaxx0 = vaxx - mean(vaxx)
		 
		 
 ## for maps 
 
   co.sp3$NAME_2[c(1119, 1125, 2671)-1]  # "Dukes"     "Nantucket" "San Juan" 
   co.sp3$index = 1:2831
   co.sp3a = subset(co.sp3, !(index %in% (c(1119, 1125, 2671)-1)) )
   co.sp3a$deth =  pmin(totDeth[-1095],500)                 # totDeth[-1095]    cnew3$dethFeb[-1095]
   spplot(co.sp3a,"deth", main="Tot deaths per 100K")
   
		# reveals some large outliers 

      ww = which(totDeth > 500)		
	  cnew3[ww,1:20]
	  
	  
	  co.sp3a$deth = pmin(totDeth[-1095], 500)
      spplot(co.sp3a,"deth")
	   
	    yhatper = yhat$x/cnew3$TotalPop*1e5        # for plotting predicted values from Binomial regression
        co.sp3a$deth =  pmin(yhatper[-1095],500)   
		spplot(co.sp3a,"deth", main="Predicted deaths per 100K")
	   
	   
	   ## for CARBayes?? 
	   
	   st = cnew3$state
       Popul.mil = cnew3$TotalPop/1e+6
	   
		Ndeath = round(totDeth* cnew3$TotalPop/1e5)    # awkward, now have to go back to counts
		
		formula <-  Ndeath ~ cnew3$vaxx0 + st + Popul.mil
     model <- S.CARleroux(formula=formula, family="binomial", trials= cnew3$TotalPop, W = co.mat3a, burnin=20000, n.sample=100000, thin = 100)
	 
	 
	  ## prep the state vars for MCMC
	  
	  st.names = levels(as.factor(cnew3$state))
	  nst = length(st.names)
	  mem = NULL 
	  for (i in 1:nst){
	    mem[[i]] = which(cnew3$state == st.names[i])	  	  
	  }
	 
	  
	   
	  
	 
    
	 
	   
	  ### =================== some checks ... 
	            	
	nnhb = apply(co.mat3,2,"sum")
    wt = which(sp.st == "Texas")
	
	 nnhb[2210:2463]
  [1] 5 6 8 4 7 7 8 6 6 6 6 8 5 7 7 7 7 7 6 5 6 4 6 7 7 5 7 6 5 6 2 5 8 5 6 4
 [37] 7 5 6 7 5 6 6 7 5 5 5 6 5 5 7 5 9 8 5 5 6 6 7 6 6 5 8 5 7 7 7 6 7 3 7 7
 [73] 5 6 8 7 7 6 5 7 5 7 6 3 7 6 6 5 8 7 7 4 6 6 7 6 6 6 6 6 7 6 6 6 5 7 8 5
[109] 6 8 5 7 7 7 5 8 6 4 7 6 7 6 5 5 6 7 7 7 6 6 4 7 6 7 7 5 4 6 6 7 7 6 6 6
[145] 7 7 7 5 7 5 5 8 7 6 5 7 6 4 5 6 6 6 6 6 6 6 5 7 6 6 8 7 7 5 5 7 7 4 7 6
[181] 6 7 5 6 5 7 6 6 3 4 7 7 4 9 6 6 8 7 4 6 7 5 5 5 6 7 7 7 6 7 7 7 4 4 5 6
[217] 7 5 6 6 6 4 7 7 4 8 6 5 4 7 6 8 5 6 6 6 6 6 7 7 6 7 6 7 3 6 5 5 6 7 5 6
[253] 3 6

  > sp.co[2213:2466]
  [1] "Anderson"(5)      "Andrews"(6)       "Angelina"(8)     "Aransas" (4)     
  [5] "Archer"        "Armstrong"     "Atascosa"      "Austin"       
  [9] "Bailey"        "Bandera"       "Bastrop"       "Baylor"       
 [13] "Bee"           "Bell"          "Bexar"         "Blanco"       
 [17] "Borden"        "Bosque"        "Bowie"         "Brazoria"     
 [21] "Brazos"        "Brewster"      "Briscoe"       "Brooks"       
 [25] "Brown"         "Burleson"      "Burnet"        "Caldwell"     
 [29] "Calhoun"       "Callahan"      "Cameron"       "Camp"         
 [33] "Carson"        "Cass"          "Castro"        "Chambers"     
 [37] "Cherokee"      "Childress"     "Clay"          "Cochran"      
 [41] "Coke"          "Coleman"       "Collin"        "Collingsworth"
 [45] "Colorado"      "Comal"         "Comanche"      "Concho"       
 [49] "Cooke"         "Coryell"       "Cottle"        "Crane"        
 [53] "Crockett"      "Crosby"        "Culberson"     "Dallam"       
 [57] "Dallas"        "Dawson"        "Deaf Smith"    "Delta"        
 [61] "Denton"        "Dewitt"        "Dickens"       "Dimmit"       
 [65] "Donley"        "Duval"         "Eastland"      "Ector"        
 [69] "Edwards"       "El Paso"       "Ellis"         "Erath"        
 [73] "Falls"         "Fannin"        "Fayette"       "Fisher"       
 [77] "Floyd"         "Foard"         "Fort Bend"     "Franklin"     
 [81] "Freestone"     "Frio"          "Gaines"        "Galveston"    
 [85] "Garza"         "Gillespie"     "Glasscock"     "Goliad"       
 [89] "Gonzales"      "Gray"          "Grayson"       "Gregg"        
 [93] "Grimes"        "Guadalupe"     "Hale"          "Hall"         
 [97] "Hamilton"      "Hansford"      "Hardeman"      "Hardin"       
[101] "Harris"        "Harrison"      "Hartley"       "Haskell"      
[105] "Hays"          "Hemphill"      "Henderson"     "Hidalgo"      
[109] "Hill"          "Hockley"       "Hood"          "Hopkins"      
[113] "Houston"       "Howard"        "Hudspeth"      "Hunt"         
[117] "Hutchinson"    "Irion"         "Jack"          "Jackson"      
[121] "Jasper"        "Jeff Davis"    "Jefferson"     "Jim Hogg"     
[125] "Jim Wells"     "Johnson"       "Jones"         "Karnes"       
[129] "Kaufman"       "Kendall"       "Kenedy"        "Kent"         
[133] "Kerr"          "Kimble"        "King"          "Kinney"       
[137] "Kleberg"       "Knox"          "La Salle"      "Lamar"        
[141] "Lamb"          "Lampasas"      "Lavaca"        "Lee"          
[145] "Leon"          "Liberty"       "Limestone"     "Lipscomb"     
[149] "Live Oak"      "Llano"         "Loving"        "Lubbock"      
[153] "Lynn"          "Madison"       "Marion"        "Martin"       
[157] "Mason"         "Matagorda"     "Maverick"      "McCulloch" 
	
	
	co.spT <- subset(co.sp3, NAME_1 =="Texas")
	co.nbT <- poly2nb(co.spT)
	 nonhbT =  which(card(co.nbT) == 0)
	   co.matT <- nb2mat(co.nbT, style="B")
	 nnhbT = apply(co.matT,2,"sum")           # checks out!
	 
  [1] 5 5 8 4 7 7 8 6 4 6 6 8 5 7 7 7 7 7 3 5 6 4 6 7 7 5 7 6 5 6 2 5 8 3 6 4
 [37] 7 4 4 5 5 6 6 5 5 5 5 6 4 5 7 5 9 8 3 3 6 6 5 6 6 5 8 5 7 7 7 6 7 1 7 7
 [73] 5 5 8 7 7 6 5 7 5 7 5 3 7 6 6 5 8 7 4 4 6 6 7 6 6 5 4 6 7 5 4 6 5 5 8 5
[109] 6 8 5 7 7 7 4 8 6 4 7 6 7 6 4 5 6 7 7 7 6 6 4 7 6 7 7 5 4 6 6 5 7 6 6 6
[145] 7 7 7 3 7 5 3 8 7 6 4 7 6 4 5 6 6 6 6 6 6 6 5 7 4 6 8 7 7 5 5 3 7 4 5 5
[181] 4 7 3 6 4 7 6 6 3 4 7 7 4 7 5 6 8 7 4 6 7 4 5 5 6 7 7 7 6 5 5 7 4 4 5 6
[217] 7 5 6 6 6 4 7 7 4 8 6 5 4 7 6 8 5 6 6 6 6 6 7 7 6 5 4 5 3 6 5 4 6 7 4 6
[253] 3 6





## this code below is obsolete, as Badr file has basically the same data as ours.
## ============================   CODE for reading the Badr file "COVID-19.csv" 

# US26147,2022-02-01,38610,0,Confirmed,Total,Total,JHU    this is the format 
# US26147,2022-02-01,38638,0,Confirmed,Total,Total,NYT
# US26147,2022-02-01,751,0,Deaths,Total,Total,JHU
# US26147,2022-02-01,754,0,Deaths,Total,Total,NYT
# US26147,2022-02-02,38890,280,Confirmed,Total,Total,JHU
# US26147,2022-02-02,38914,276,Confirmed,Total,Total,NYT


# Open the file for reading
file_path = "C:/localH/research/Awal/COVID-19.csv"
fin = file(file_path, "r")
fout <- file("C:/localH/research/Awal/COVID-19us.csv", "w")
nlines = 0   # the file has  30,963,308 lines!

# Read and process lines one by one
 while(TRUE) {        # while(nlines < 100) {
  linea <- readLines(fin, n = 1)
  if (length(linea) == 0) {
    break  # Exit loop if end of file is reached
  } 
  nlines = nlines + 1
  lxx = strsplit(linea,',')[[1]]
  str1 = lxx[1]
  str2 = lxx[2]
  str5 = lxx[5]
  str8 = lxx[8]
  
  if (startsWith(str1, "US") & startsWith(str2, "2021") &(str5 == "Deaths") & (str8 == "JHU") ) {
    # Write the line to the output file
    cat(linea, "\n", file = fout)
  }
  
 }
  
  # Close the file connection
close(fin)
close(fout)
  
  #================================
  ## the US file is still too big, so here is a refinement 
  
  file_path = "C:/localH/research/Awal/COVID-19us.csv"
fin = file(file_path, "r")
fout <- file("C:/localH/research/Awal/COVID-19usDeaths2021NYT.csv", "w")
nlines = 0

# Read and process lines one by one
 while(TRUE) {        # while(nlines < 100) {
  linea <- readLines(fin, n = 1)
  if (length(linea) == 0) {
    break  # Exit loop if end of file is reached
  } 
  nlines = nlines + 1
  lxx = strsplit(linea,',')[[1]]
  str1 = lxx[1]
  str2 = lxx[2]
  str5 = lxx[5]
  str8 = lxx[8]
  
  if (startsWith(str2, "2021") &(str5 == "Deaths") & startsWith(str8, "NYT")) {
    # Write the line to the output file
    cat(linea, "\n", file = fout)
  }
  
 }
  
  # Close the file connection
close(fin)
close(fout)

 cntyUS = read.csv("C:/localH/research/Awal/COVID-19usDeaths2021NYT.csv", header=F)  
 cntyUS = cntyUS[,1:4]
  nms = levels(as.factor(cntyUS[,1]))
  
  cid = "US31111"   # Miami-Dade "US01095"  
  subs = (cntyUS[,1] == cid)
  Yd = cntyUS[subs,3]
  plot(Yd,type="l")
  
  
  
  
  
	
	
	
	

    
   
  



