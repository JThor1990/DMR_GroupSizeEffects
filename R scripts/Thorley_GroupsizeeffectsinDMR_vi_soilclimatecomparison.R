#----------------------------------------------------------------
# Damaraland mole-rat breeders do not rely on helpers for reproduction or survival
#
# vi- "Comparison of the soil properties and climate characteristics across three medium to long-term DMR study populations"
#
# R script
# Authors: Jack Thorley, Hanna Bensch
# Contact: jack.thorley1@gmail.com; jbt27@cam.ac.uk
#----------------------------------------------------------------

# In the script I will compare temperature, rainfall, and soil properties at three study sites: KRR (our study site), Tswalu (Southern Kalahari, South Africa), and Dordabis (Namibia). 

# Rainfall and temperature data come from https://power.larc.nasa.gov/data-access-viewer/, and cover the period from 1st Jan 1981 to 31st Dec 2020

# Soil data are extracted from the soilGrids database, accessible through R via the 'soilDB' package (see below)

#---------------------------

library(tidyverse) ; library(ggmap) ; library(raster); library(patchwork) ; library(ggrepel) ; library(soilDB) ;  library(multcomp)


# store the coordinates for later use
KRR_lat <- -26.97856 ; KRR_long <- 21.79000
Tswalu_lat <- -27.43337 ; Tswalu_long <- 22.26670
Dordabis_lat <- -22.96667 ; Dordabis_long <- 17.68364

# read in climate data 
climate <- read.csv("FieldMR_ClimateSoil.csv", header = TRUE) 

#----------------------------------
#  Making a map for the DMR distribution

loc <- getwd() # All the DMR shape files and components need to be in current working directory ("DMR shapefiles")
DMR_dist <- paste0(loc, "data_0.shp")
DMR_dist <- shapefile(DMR_dist)
DMR_dist <- DMR_dist[1,]
  
kalahari <- paste0(loc, "kalaharidesert_gai.shp")
kalahari <- shapefile(kalahari)

#Set your API Key
register_google(key = "########") # Need to set up own API here
has_google_key() # Check = TRUE

p <- ggmap(get_googlemap(center = c(lon = 21.832459, lat = -23.68560),                  
                         maptype = "satellite", zoom = 6))

world <- map_data("world") %>%  # we already did this, but we can do it again
  filter(region != "Lesotho")

map1 <- p +   # N.B. warning relates to the eclusion of the rest of the world map country boundaries
  geom_path(data = world, aes(x = long, y = lat, group  = group), colour = "black") + 
  geom_polygon(data = kalahari, aes(x = long, y = lat), alpha = 0.5, fill = "darkorange", linewidth = 1) + 
  geom_polygon(data = DMR_dist, aes(x = long, y = lat), alpha = 0,  colour = "white", lty = 2, lwd = 1.2) + 
  xlab('Longitude') + 
  ylab('Latitude') +
  theme_bw() + 
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  geom_point(aes(x = KRR_long, y = KRR_lat), pch = 24, fill = "red", size = 3) + 
  geom_point(aes(x = Dordabis_long, y = Dordabis_lat), pch = 21, fill = "red", size = 3) +
  geom_point(aes(x = Tswalu_long, y = Tswalu_lat), pch = 21, fill = "red", size = 3) + 
  geom_text(aes(x = KRR_long + 0.6, y = KRR_lat + 0.3), label = "KRR") + 
  geom_text(aes(x = Tswalu_long + 0.9, y = Tswalu_lat), label = "Tswalu") + 
  geom_text(aes(x = Dordabis_long + 0.6, y = Dordabis_lat + 0.3), label = "Dordabis") + 
  geom_text(aes(x = 22.5, y = -23.5), label = "KALAHARI \n DESERT", colour = "black") + 
  geom_text(aes(x = 25, y = -28.5), label = "SOUTH AFRICA", colour = "white", size = 5) +
  geom_text(aes(x = 24, y = -19.5), label = "BOTSWANA", colour = "white", size = 5) + 
  geom_text(aes(x = 17, y = -25), label = "NAMIBIA", colour = "white", size = 5)
  

#--------------------------------------------------

# CLIMATE:

# Plot the average monthly climate at each site

#Points represent the mean monthly maximum temperature (black), minimum temperature (grey) and rainfall across the period, with between-year error bars noted by the standard error bars. Points are scaled to the median within-month variance of each climatic variable. 

# tidy the data to get all the information
sem <-  function(x) sd(x)/sqrt(length(x))

climatesummary <- climate %>% 
  group_by(SITE, MO, YEAR) %>%  # get the average mean monthly temperature
  summarise(totalrainfall = sum(PRECIPTATION_mmday), 
            meantemp = mean(T2M), maxtemp = max(T2M_MAX), 
            mintemp = min(T2M_MIN)) %>% 
  group_by(SITE, MO) %>% 
  summarise(meanmonthlyrain = mean(totalrainfall), semmonthlyrain = sem(totalrainfall),
            meanavgtemp = mean(meantemp), semavgtemp = sem(meantemp), 
            meanmaxtemp  = mean(maxtemp), semmaxtemp = sem(maxtemp), 
            meanmintemp  = mean(mintemp), semmintemp = sem(mintemp)) %>% 
  mutate(month = month.abb, 
         rainupper = meanmonthlyrain + semmonthlyrain, 
         rainlower = meanmonthlyrain - semmonthlyrain, 
         meantempupper = meanavgtemp + semavgtemp,
         meantemplower = meanavgtemp - semavgtemp,
         maxtempupper = meanmaxtemp + semmaxtemp,
         maxtemplower = meanmaxtemp - semmaxtemp,
         mintempupper = meanmintemp + semmintemp,
         mintemplower = meanmintemp - semmintemp) %>% 
data.frame()
climatesummary$month <- factor(climatesummary$month, levels = month.abb)

# generate the rainfall plot for all sites 
allsites_rainplot <- ggplot(data = climatesummary, 
                            aes(x = month, y = meanmonthlyrain, group = SITE)) + 
  geom_path(colour = "lightgrey") +
  geom_linerange(aes(ymin = rainlower , ymax = rainupper), colour = "blue") + 
  geom_point(aes(size = semmonthlyrain), colour = "blue") +
  facet_wrap(~SITE) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 13),
        legend.position = "none", 
        strip.background = element_rect(fill = "goldenrod"), 
        strip.text = element_text(size = 14), 
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(y = expression(paste("Mean monthly rainfall (mm ", phantom(.) %+-% phantom(.),"SEM)")))


# generate the temperature plot for all sites
allsites_tempplot <- ggplot(data = climatesummary) + 
  #mean temps
  geom_line(aes(x = month, y = meanavgtemp, group = 1), colour = "black") + 
  geom_linerange(aes(x = month, y = meanavgtemp, ymax = meantempupper, ymin = meantemplower), colour = "black") + 
  geom_point(aes(x = month, y = meanavgtemp, size = semavgtemp), colour = "black") +
  #max temps
  geom_line(aes(x = month, y = meanmaxtemp, group = 1), colour = "firebrick") + 
  geom_linerange(aes(x = month, y = meanmaxtemp, ymax = maxtempupper, ymin = maxtemplower), colour = "firebrick") + 
  geom_point(aes(x = month, y = meanmaxtemp, size = semmaxtemp), colour = "firebrick") +
  #min temps
  geom_line(aes(x = month, y = meanmintemp, group = 1), colour = "dodgerblue") + 
  geom_linerange(aes(x = month, y = meanmintemp, ymax = mintempupper, ymin = mintemplower), colour = "dodgerblue") + 
  geom_point(aes(x = month, y = meanmintemp, size = semmintemp), colour = "dodgerblue") +
  facet_wrap(~SITE) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14),
        legend.position = "none", 
        strip.background = element_rect(fill = "goldenrod"), 
        strip.text = element_text(size = 14), 
        panel.grid = element_blank()) +
  labs(x = "Month", y = "Temperature (°C)") + 
  scale_size(range = c(1,5))   # note that the temperature s.e.m. are very small so cannot be seen on the plot

climateplots <- allsites_rainplot / allsites_tempplot   # makes use of patchwork syntax


# Plot the average annual climate at each site

# Total annual rainfall (mm)
annualrain <- climate %>% 
  group_by(SITE, YEAR) %>% 
  summarise(totalannualrainfall = sum(PRECIPTATION_mmday)) %>% 
  filter(YEAR > 1981) %>% 
  data.frame()

quantiles <- annualrain %>% 
  group_by(SITE) %>% 
  summarise(mean = mean(totalannualrainfall), 
            l95 = quantile(totalannualrainfall, probs = c(0.025)), 
            u95 = quantile(totalannualrainfall, probs = c(0.975)))
  
climateplots2 <- ggplot(annualrain, aes(x = totalannualrainfall, group = SITE, fill = SITE, colour = SITE)) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 15) + 
  geom_vline(data = quantiles, aes(xintercept = mean, colour = SITE), linetype = 1, linewidth = 1.5) +
  geom_vline(data = quantiles, aes(xintercept = l95, colour = SITE), linetype = 2, linewidth = 1.2) +
  geom_vline(data = quantiles, aes(xintercept = u95, colour = SITE), linetype = 2, linewidth = 1.2) +
  facet_wrap(~SITE, ncol = 1) + 
  theme_bw() +  
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "goldenrod"), 
        strip.text = element_text(size = 12), 
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Count", x = "Total annual rainfall (mm)") + 
  scale_x_continuous(breaks = seq(100, 1800, 100))

# Annual temperature (mean daily maximum temperature)
annualtemp <- climate %>% 
  group_by(SITE, YEAR) %>% 
  summarise(meandailytemperature = mean(T2M_MAX)) %>% 
  filter(YEAR > 1981) %>% 
  data.frame()

quantiles2 <- annualtemp %>% 
  group_by(SITE) %>% 
  summarise(mean = mean(meandailytemperature), 
            l95 = quantile(meandailytemperature, probs = c(0.025)), 
            u95 = quantile(meandailytemperature, probs = c(0.975)))

climateplots3 <- ggplot(annualtemp, aes(x = meandailytemperature, group = SITE, fill = SITE, colour = SITE)) + 
  geom_histogram(position = 'identity', alpha = 0.4, bins = 15) + 
  geom_vline(data = quantiles2, aes(xintercept = mean, colour = SITE), linetype = 1, linewidth = 1.5) +
  geom_vline(data = quantiles2, aes(xintercept = l95, colour = SITE), linetype = 2, linewidth = 1.2) +
  geom_vline(data = quantiles2, aes(xintercept = u95, colour = SITE), linetype = 2, linewidth = 1.2) +
  facet_wrap(~SITE, ncol = 1) + 
  theme_bw() +  
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "goldenrod"), 
        strip.text = element_text(size = 12), 
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        panel.grid = element_blank()) +
  labs(y = "Count", x = "Mean daily maximum \ntemperature (°C)")

climateplots2 + 
  labs(tag = "(a)", size = 14) + 
  climateplots3 + 
  labs(tag = "(b)", size = 14) 


# compare the means statistically
annualrain$SITE <- factor(annualrain$SITE, levels = c("KRR", "Dordabis", "Tswalu"))

m1 <- lm(totalannualrainfall ~ SITE, data = annualrain)
summary(m1)
anova(m1)

cont1 <- glht(m1, linfct = mcp(SITE = "Tukey")) # get the contrasts among the different levels for site
summary(cont1)

annualtemp$SITE <- factor(annualtemp$SITE, levels = c("KRR", "Dordabis", "Tswalu"))
m2 <- lm(meandailytemperature ~ SITE, data = annualtemp)
summary(m2)
anova(m2)

cont2 <- glht(m2, linfct = mcp(SITE = "Tukey"))
summary(cont2)

#--------------------------------

# SOIL PROPERTIES:


# To get rid of any of uncertainty in the soilGrid estimated properties at a single grid cell, 
# I will create a 2km by 2km bounding box around each location and sample the properties at these locations too. 
# (Note that SoilGrid provides data at a 250m by 250m resolution, so the points will be distinct)
# This will generate four additional data points for each site. I will then take the average soil properties across each of the five points as representative of the general soil properties of each site. 

soil_latlon <- data.frame(id = c("KRR", "Tswalu", "Dordabis"), 
                          lat = c(KRR_lat, Tswalu_lat, Dordabis_long), 
                          lon = c(KRR_long, Tswalu_long, Dordabis_long))

# need to 

coordinates(soil_latlon) <- c("lon", "lat")
proj4string(soil_latlon) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
soil_latlon <- spTransform(soil_latlon, CRS("+proj=utm +zone=34 +south +ellps=WGS84"))
soil_latlon <- data.frame(soil_latlon) %>%  dplyr::select(-optional)

# expand to get the coordinates in a the corners of a 2km bounding box
fourcorner <- function(x, y) {
  
  p1_x <- x - 2000/2 # 1000m to the left
  p1_y <- y - 2000/2 # 1000m down

  p2_x <- x - 2000/2 # etc...
  p2_y <- y + 2000/2
  
  p3_x <- x + 2000/2
  p3_y <- y + 2000/2
  
  p4_x <- x + 2000/2
  p4_y <- y - 2000/2
  
 return(data.frame(lon = c(p1_x, p2_x, p3_x, p4_x), 
                   lat = c(p1_y, p2_y, p3_y, p4_y))) 
  
}

fourcorner(soil_latlon$lon[1], soil_latlon$lat[1]) # shows function works. 

soil_latlon2 <- fourcorner(soil_latlon$lon, soil_latlon$lat) %>% 
  mutate(id = rep(c("KRR", "Tswalu", "Dordabis"), times = 4)) %>% 
  dplyr::select(id, lon, lat)

soil_latlon <- rbind(soil_latlon, soil_latlon2) # back to lat lon rather than UTM for the data extraction
coordinates(soil_latlon) <- c("lon", "lat")
proj4string(soil_latlon) <- CRS("+proj=utm +zone=34 +south +ellps=WGS84")
soil_latlon <- spTransform(soil_latlon, CRS("+proj=longlat +datum=WGS84 +no_defs"))
soil_latlon <- data.frame(soil_latlon) %>%  dplyr::select(-optional, id, lat, lon)

soildf <- fetchSoilGrids(   # extract the soil properties for each coordinate
  soil_latlon %>% 
    mutate(id = 1:15),  # needs a unique id to work
  loc.names = c("id", "lat", "lon"),
  verbose = TRUE,
  progress = TRUE)

soilfinal <- data.frame(label = soildf$label, 
           id = soil_latlon$id, 
           bdodmean = soildf$bdodmean, # bulk density of fine earth fraction
           cfvomean = soildf$cfvomean, # volumetric fraction of coarse fragments
           claymean = soildf$claymean, # proportion of clay particles
           nitrogenmean = soildf$nitrogenmean, # total Nitrogen
           sandmean = soildf$sandmean, # proportion of sand particles
           siltmean = soildf$siltmean, # proportion of silt particles
           socmean = soildf$socmean,  #soil organic carbon
           soildf$hzdept, # horizon depth upper (not used)
           soildf$hzdepb) # horizon depth lower (not used)

# Tidy the data to get a summary
soilfinal <- soilfinal %>% 
  filter(label != "100-200") %>% # only get the upper 1m of soil, the typical depth of mole-rat burrows (though including makes little quantitative difference)
  group_by(id) %>% 
  summarise(bdod = mean(bdodmean),
            cfvoc = mean(cfvomean),
            clay = mean(claymean),
            nitrogenmean = mean(nitrogenmean), 
            sandmean = mean(sandmean), 
            siltmean = mean(siltmean), 
            socmean = mean(socmean)) %>% 
  data.frame()
            
# Comparison across the species' range

# Lastly, I want to compare whether the soil properties are representative of the entire DMR distribution. 
# For this, I will take regular points from the DMR distribution at 50km intervals
DMR_dist_UTM <- spTransform(DMR_dist, CRS("+proj=utm +zone=34 +south +ellps=WGS84"))
soil_samples <- spsample(DMR_dist_UTM, type = "regular", cellsize = 50000) # 50000m intervals
soil_samples <- spTransform(soil_samples, CRS("+proj=longlat +datum=WGS84 +no_defs"))
pts_to_plot <- data.frame(soil_samples)
names(pts_to_plot) <- c("lon", "lat")  # 350 samples across the range

map2 <- p +   
  geom_path(data = world, aes(x = long, y = lat, group  = group), colour = "black") + 
  geom_polygon(data = kalahari, aes(x = long, y = lat), alpha = 0.5, fill = "darkorange", linewidth = 1) + 
  geom_polygon(data = DMR_dist, aes(x = long, y = lat), alpha = 0,  colour = "white", lty = 2, lwd = 1.2) + 
  geom_point(data = pts_to_plot, aes(x = lon, y = lat),  colour = "cyan", pch = 16) + 
  geom_point(data = pts_to_plot, aes(x = lon, y = lat),  colour = "black", pch = 1) +
  xlab('Longitude') + 
  ylab('Latitude') +
  theme_bw() + 
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  geom_point(aes(x = KRR_long, y = KRR_lat), pch = 24, fill = "red", size = 3) + 
  geom_point(aes(x = Dordabis_long, y = Dordabis_lat), pch = 21, fill = "red", size = 3) +
  geom_point(aes(x = Tswalu_long, y = Tswalu_lat), pch = 21, fill = "red", size = 3) + 
  geom_text(aes(x = KRR_long + 0.6, y = KRR_lat + 0.3), label = "KRR") + 
  geom_text(aes(x = Tswalu_long + 0.9, y = Tswalu_lat), label = "Tswalu") + 
  geom_text(aes(x = Dordabis_long + 0.6, y = Dordabis_lat + 0.3), label = "Dordabis") + 
  geom_text(aes(x = 22.5, y = -23.5), label = "KALAHARI \n DESERT", colour = "black") + 
  geom_text(aes(x = 25, y = -28.5), label = "SOUTH AFRICA", colour = "white", size = 5) +
  geom_text(aes(x = 24, y = -19.5), label = "BOTSWANA", colour = "white", size = 5) + 
  geom_text(aes(x = 17, y = -25), label = "NAMIBIA", colour = "white", size = 5)

soildf2 <-  data.frame(soil_samples) %>% 
  mutate(id = 1:n()) %>% 
  rename(lon = x1, lat = x2) %>% 
  dplyr::select(id, lat, lon)
  
soildf_random <- fetchSoilGrids(
  soildf2,
  loc.names = c("id", "lat", "lon"),
  verbose = TRUE,
  progress = TRUE) # will take some time

# collate into data.frame
soilfinal2 <- data.frame(label = soildf_random$label, 
                        id = soildf_random$id, 
                        bdodmean = soildf_random$bdodmean, # bulk density of fine earth fraction
                        cfvomean = soildf_random$cfvomean, # volumetric fraction of coarse fragments
                        claymean = soildf_random$claymean, # proportion of clay particles
                        nitrogenmean = soildf_random$nitrogenmean, # total Nitrogen
                        sandmean = soildf_random$sandmean, # proportion of sand particles
                        siltmean = soildf_random$siltmean, # proportion of silt particles
                        socmean = soildf_random$socmean,  #soil organic carbon
                        hzdpt = soildf_random$hzdept, # horizon depth upper (not used)
                        hzdepb = soildf_random$hzdepb) # horizon depth lower (not used)

soilfinal2 <- soilfinal2 %>% 
  filter(label != "100-200") %>% # filter to less than 1m
  group_by(id) %>% 
  summarise(bdod = mean(bdodmean, na.rm = TRUE),
            cfvoc = mean(cfvomean, na.rm = TRUE),
            clay = mean(claymean, na.rm = TRUE),
            nitrogenmean = mean(nitrogenmean, na.rm = TRUE), 
            sandmean = mean(sandmean, na.rm = TRUE), 
            siltmean = mean(siltmean, na.rm = TRUE), 
            socmean = mean(socmean, na.rm = TRUE)) %>% 
  data.frame()

# Plot the variation in soil properties as a PCA
pca_data <- data.frame(rbind(soilfinal, soilfinal2))
row.names(pca_data) <- pca_data$id
pca_data <- na.omit(pca_data)  # 7 NAs- which we can remove for the purpose of the visualisation

pc <- prcomp(x = pca_data[-1],
             center = TRUE, 
             scale. = TRUE)
print(pc)
summary(pc)

pca_plotdf <- data.frame(id = row.names(pc$x), pc$x)

p3 <- ggplot(pca_plotdf, aes(x = PC1, y = PC2)) + 
  geom_hline(yintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, colour = "darkgrey") +
  geom_point(pch = 21, col = "black", fill = "cyan") +
  theme_bw() + 
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  labs(x = "PC1 (42.8% var)", y = "PC2 (20.0% var)")

datapc <- data.frame(varnames = rownames(pc$rotation), pc$rotation)
datapc$varnames <- c("Bulk\ndensity", "Volumetric\nfraction", "% Clay", "Nitrogen", "% Sand", "% Silt", "Organic\ncarbon")
mult <- min(
  (max(pca_plotdf[,"PC2"]) - min(pca_plotdf[,"PC2"])/(max(datapc[,"PC2"]) - min(datapc[,"PC2"]))),
  (max(pca_plotdf[,"PC1"]) - min(pca_plotdf[,"PC1"])/(max(datapc[,"PC1"]) - min(datapc[,"PC1"])))
)

datapc <- transform(datapc,
                    v1 = 0.7 * mult * PC1,
                    v2 = 0.7 * mult * PC2) # transform the aid the visualisation

datapc_text <- datapc 
datapc_text$v1[c(3,6)] <- datapc_text$v1[c(3,6)] - 0.9
datapc_text$v1[c(1, 5)] <- datapc_text$v1[c(1, 5)] + 0.9
datapc_text$v2[c(1,2)] <- datapc_text$v2[c(1,2)] - 0.9
datapc_text$v2[c(4,7)] <- datapc_text$v2[c(4,7)] + 0.75
datapc_text$v2[6] <- datapc_text$v2[6] - 0.3
  
p3 <- p3 + 
  geom_segment(data = datapc, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.2,"cm")), 
               alpha = 0.75, linewidth = 1.2, color = "purple") + 
  geom_label(data = datapc_text, aes(x = v1, y = v2, label = varnames), size = 4, color = "purple") +
  geom_point(data = filter(pca_plotdf, id == "KRR"), pch = 24, fill = "red", size = 3) + 
  geom_point(data = filter(pca_plotdf, id == "Tswalu"), pch = 21, fill = "red", size = 3) +
  geom_point(data = filter(pca_plotdf, id == "Dordabis"), pch = 21, fill = "red", size = 3) +
  geom_label_repel(data = filter(pca_plotdf, id %in% c("KRR", "Tswalu", "Dordabis")), 
                  force = 40, aes(label = id), size = 4) + 
  xlim(c(-6, 6)) + 
  ylim(c(-6, 4)) 


#----------------------------------

# join the two plots together

plot_soil <- (map2 + 
  labs(tag = "(a)") + 
  theme(plot.tag = element_text(size = 14)))  +
  (p3 + 
  labs(tag = "(b)") + 
  theme(plot.tag = element_text(size = 14)))
  

#####-----------------  END   ------------------####
