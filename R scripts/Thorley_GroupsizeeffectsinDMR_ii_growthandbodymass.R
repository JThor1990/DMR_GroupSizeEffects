#----------------------------------------------------------------
# Damaraland mole-rat breeders do not rely on helpers for reproduction or survival
#
# ii- "Early life growth and adult body mass" 
#
# R script
# Authors: Jack Thorley, Hanna Bensch
# Contact: jack.thorley1@gmail.com; jbt27@cam.ac.uk
#----------------------------------------------------------------

library(tidyverse) ; library(nlme) ; library(MASS) ; library(multcomp) 

# In this script we model growth in body mass and skeletal size (via teeth width), and examine adult body mass. 

# read in the appropriate data set (will need to set your own working directory here)
bodymass <- read.csv("FieldMR_Growth_BodyMass.csv", header = TRUE)   
teethwidth <- read.csv("FieldMR_Growth_Teethwidth.csv", header = TRUE) 

# separate by sex (body mass)
malebodymass <- filter(bodymass, Sex == "Male")
femalebodymass <- filter(bodymass, Sex == "Female")

# separate by sex (teeth width)
maleteeth <- filter(teethwidth, Sex == "Male")
femaleteeth <- filter(teethwidth, Sex == "Female")

# information on body mass
dim(malebodymass) # 456
dim(femalebodymass) # 381

length(unique(malebodymass$AnimalID)) # 214
length(unique(femalebodymass$AnimalID)) # 193

length(unique(malebodymass$GroupID)) # 39
length(unique(femalebodymass$GroupID)) # 48

malebodymass %>% 
  group_by(AnimalID) %>% 
  summarise(count = n()) %>% 
  summarise(mean(count), sd(count))

femalebodymass %>% 
  group_by(AnimalID) %>% 
  summarise(count = n()) %>% 
  summarise(mean(count), sd(count))

mean(femalebodymass$GroupSize) ; sd(femalebodymass$GroupSize)
mean(malebodymass$GroupSize) ; sd(malebodymass$GroupSize)

# information on teeth width
dim(maleteeth) # 381
dim(femaleteeth) # 328

length(unique(maleteeth$AnimalID)) # 198
length(unique(femaleteeth$AnimalID)) # 180

length(unique(maleteeth$GroupID)) # 38
length(unique(femaleteeth$GroupID)) # 47

maleteeth %>% 
  group_by(AnimalID) %>% 
  summarise(count = n()) %>% 
  summarise(mean(count), sd(count))

femaleteeth %>% 
  group_by(AnimalID) %>% 
  summarise(count = n()) %>% 
  summarise(mean(count), sd(count))

mean(femaleteeth$GroupSize) ; sd(femaleteeth$GroupSize)
mean(maleteeth$GroupSize) ; sd(maleteeth$GroupSize)


## MODEL GROWTH #-------------------------------------------------------------------------------------------------------
 # Using an interval equation for the von bertalanffy function from Schoener and Shoener 1978 Copeia

# BODY MASS 
#---------------

# first run a basic nls to get an idea of starting parameters
# Females body mass - baseline model (without group size)

summary(nls(weight2 ~ A - (A - weight)*exp(-k*timediff),
            data= femalebodymass,
            start=c(A=115,k=0.004))) # A = 132.2, k = 0.00263

vonbert.female <- nlme(weight2 ~ A - (A - weight)*exp(-k*timediff),
                       fixed=A+k~1, 
                       random= list(AnimalID = pdDiag(A + k ~ 1)),
                       data= femalebodymass,
                       start=c(A=115,k=0.00415), 
                       na.action = na.omit, 
                       control = nlmeControl(maxIter = 2000, pnlsTol = 0.01, tolerance = 1e-3))
summary(vonbert.female)   # A = 118.85480 ; k = 0.00415     

# males body mass - baseline model (without group size)

summary(nls(weight2 ~ A - (A - weight)*exp(-k*timediff),
            data= malebodymass,
            start=c(A=170,k=0.004))) # A = 168, k = 0.002375

vonbert.male <- nlme(weight2 ~ A - (A - weight)*exp(-k*timediff),
                     fixed=A+k~1, 
                     random= list(AnimalID = pdDiag(A + k ~ 1)),
                     data= malebodymass,
                     start=c(A=170,k=0.01), 
                     na.action = na.omit)
summary(vonbert.male)   # A = 149.9768; k = 0.0035  


# Now include group size and rainfall into the models
# First standardise group size and rainfall
mean(femalebodymass$GroupSize) ; sd(femalebodymass$GroupSize) # 10.88 # 6.47
mean(malebodymass$GroupSize) ; sd(malebodymass$GroupSize) # 12.57 # 6.01

femalebodymass$GroupSize.s <- as.numeric(scale(femalebodymass$GroupSize))
malebodymass$GroupSize.s <- as.numeric(scale(malebodymass$GroupSize))

# Female model 2 (with group size and rainfall
vonbert.female2 <- nlme(weight2 ~ (A + AGS*GroupSize.s) - ((A + AGS*GroupSize.s) - weight)*exp(-(k + kGS*GroupSize.s)*timediff),
                        fixed=A+k + AGS + kGS ~1, 
                        random= list(AnimalID = pdDiag(A + k ~ 1)),
                        data= femalebodymass,
                        start=c(A=118,k=0.02, AGS = 0, kGS = 0), 
                        na.action= na.omit)
summary(vonbert.female2)   

# Male model 2 (with group size and rainfall
vonbert.male2 <- nlme(weight2 ~ (A + AGS*GroupSize.s) - ((A + AGS*GroupSize.s) - weight)*exp(-(k + kGS*GroupSize.s )*timediff),
                          fixed=A+k + AGS + kGS ~1, 
                          random= list(AnimalID = pdDiag(A + k ~ 1)),
                          data= malebodymass,
                          start=c(A=150,k=0.02, AGS = -10, kGS = 0), 
                          na.action= na.omit)
summary(vonbert.male2)   

AIC(vonbert.female, vonbert.female2) %>% 
  mutate(AIC - min(AIC))
AIC(vonbert.male, vonbert.male2) %>% 
  mutate(AIC - min(AIC))  

summary(vonbert.female2)
summary(vonbert.male2)

# Plot rate of change in body mass growth (with group size accounted for)
#--------------------------------------------------

 # male body mass- rate of growth
summary(vonbert.male2)

vonbertratefunc <- function(x, A, k, AGs, kGs, GS) {
  
  time2 = (A + AGs*GS) - ((A + AGs*GS) - x)*exp(-(k + kGs*GS)*180)
  return(time2)
  }

# pick a colour blind friendly palette
colpal <- palette.colors(palette = "Okabe-Ito")[c(3,4,6)]

par(mfrow=c(1,2))
  plot(rateofbodymassgrowth ~ weight, data = malebodymass, pch = 1, col = adjustcolor("grey20", alpha.f = 0.6), xlab = "Body Mass (g)", ylab = "Rate of change in body mass (g/day)", bty = "l", las = 1, xlim = c(10, 200), ylim = c(-0.1, 0.4))
  axis(1, at = seq(0, 200, 25), labels = NA) ; axis(2, at = seq(-0.05, 0.35, 0.05), labels = NA)
  abline(h = 0, lty =2, col = adjustcolor("grey", alpha.f = 0.6))
  abline(v = 149.69626, lty = 1, col = adjustcolor("grey", alpha.f = 0.9))
  
  # predict rate of change at 4, 12, and 20 individuals
  small <- (4 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize)
  malesmall <- data.frame(time1 = seq(10, 161, 0.5),
                          time2 = vonbertratefunc(x = seq(10, 161, 0.5), 
                                                   A = 149.69626, 
                                                   AGs = -8.14804, 
                                                   k = 0.00368, 
                                                   kGs = 0.00084, 
                                                   GS = small)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malesmall, type = "l", lwd = 4, lty = 3, col = colpal[1])
  
  medium <- (12 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize)
  malemedium <- data.frame(time1 = seq(10, 149, 0.5),
                           time2 = vonbertratefunc(x = seq(10, 149, 0.5), 
                                                   A = 149.69626, 
                                                   AGs = -8.14804, 
                                                   k = 0.00368, 
                                                   kGs = 0.00084, 
                                                   GS = medium)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malemedium, type = "l", lwd = 4, lty = 2, col = colpal[2])
  
  
  large <- (20 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize)
  malelarge <- data.frame(time1 = seq(10, 141, 0.5),
                          time2 = vonbertratefunc(x = seq(10, 141, 0.5), 
                                                  A = 149.69626, 
                                                  AGs = -8.14804, 
                                                  k = 0.00368, 
                                                  kGs = 0.00084, 
                                                  GS = large)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malelarge, type = "l", lwd = 3, lty = 1, col = colpal[3])
  text(30, -0.05, "Males", cex = 1.2)


  # female body mass- rate of growth
  summary(vonbert.female2)
  plot(rateofbodymassgrowth ~ weight, data = femalebodymass, pch = 4, col = adjustcolor("grey20", alpha.f = 0.6), xlab = "Body Mass (g)", ylab = "Rate of change in body mass (g/day)", bty = "l", las = 1, xlim = c(10, 200), ylim = c(-0.1, 0.4))
  axis(1, at = seq(0, 200, 25), labels = NA) ; axis(2, at = seq(-0.05, 0.35, 0.05), labels = NA)
  abline(h = 0, lty =2, col = adjustcolor("grey", alpha.f = 0.6))
  abline(v = 118.26449, lty = 1, col = adjustcolor("grey", alpha.f = 0.9))
  
  # at 4, 12, and 20 individuals
  small <- (4 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize)
  femalesmall <- data.frame(time1 = seq(10, 129, 0.5),
                            time2 = vonbertratefunc(x = seq(10, 129, 0.5), 
                                                     A = 118.26449, 
                                                     AGs = -9.41884, 
                                                     k = 0.00447, 
                                                     kGs = 0.00155, 
                                                     GS = small)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalesmall, type = "l", lwd = 3, lty = 3, col = colpal[1])
  
  medium <- (12 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize)
  femalemedium <- data.frame(time1 = seq(10, 117, 0.5),
                             time2 = vonbertratefunc(x = seq(10, 117, 0.5), 
                                                     A = 118.26449, 
                                                     AGs = -9.41884, 
                                                     k = 0.00447, 
                                                     kGs = 0.00155, 
                                                     GS = medium)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalemedium, type = "l", lwd = 3, lty = 2, col = colpal[2])
  
  
  large <- (20 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize)
  femalelarge <- data.frame(time1 = seq(10, 105, 0.5),
                            time2 = vonbertratefunc(x = seq(10, 105, 0.5), 
                                                    A = 118.26449, 
                                                    AGs = -9.41884, 
                                                    k = 0.00447, 
                                                    kGs = 0.00155, 
                                                    GS = large)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalelarge, type = "l", lwd = 3, lty = 1, col = colpal[3])
  text(30, -0.05, "Females", cex = 1.1)


# Plot age-related variation in body mass growth by deriving the rate of change from early life 
# pup mass at birth is on average 10g, so we can forecast the rate of change in mass from this point at different group sizes
  
# get the male prediction first (at group sizes of 4, 12, and 20)
timevals <- seq(0, 730, 1)

pframe.malemass <- expand.grid(weight = 10, timediff = timevals, 
                       GroupSize.s = c((4 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize),
                                       (12 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize), 
                                       (20 - mean(malebodymass$GroupSize))/sd(malebodymass$GroupSize)))
pframe.femalemass <- expand.grid(weight = 10, timediff = timevals, 
                       GroupSize.s = c((4 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize),
                                       (12 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize), 
                                       (20 - mean(femalebodymass$GroupSize))/sd(femalebodymass$GroupSize)))
pframe.malemass$weight2 <- as.vector(predict(vonbert.male2, newdata = pframe.malemass, level = 0))
pframe.femalemass$weight2 <- as.vector(predict(vonbert.female2, newdata = pframe.femalemass, level = 0))

# now get the bootstrapped values for the mean
nresamp <- 1000
## pick new parameter values by sampling from multivariate normal distribution based on fit
pars.picked1 <- mvrnorm(nresamp, mu = fixef(vonbert.male2), Sigma = vcov(vonbert.male2))
pars.picked2 <- mvrnorm(nresamp, mu = fixef(vonbert.female2), Sigma = vcov(vonbert.female2))

## utility function
get_CI <- function(y,pref="") {
  r1 <- t(apply(y,1,quantile,c(0.025,0.975)))
  setNames(as.data.frame(r1),paste0(pref,c("lwr","upr")))
}

set.seed(101)
yvals1 <- matrix(nrow = length(timevals)*3, ncol = 1000)
yvals2 <- yvals1

vonbertratefunc2 <- function(A, k, AGS, kGS, GroupSize.s, timediff) {
  time2 = (A + AGS*GroupSize.s) - ((A + AGS*GroupSize.s) - 10)*exp(-(k + kGS*GroupSize.s)*timediff)
  return(time2)
}

# male prediction
for (i in 1:length(pframe.malemass$timediff))
{
  yvals1[i,] <- vonbertratefunc2(timediff = pframe.malemass$timediff[i], 
                                  GroupSize.s = pframe.malemass$GroupSize.s[i],
                                  A = pars.picked1[,1], k = pars.picked1[,2], 
                                  AGS = pars.picked1[,3], kGS = pars.picked1[,4])
}
ci.males <- get_CI(yvals1)
pframe.malemass <- cbind(pframe.malemass, ci.males)  

# female prediction
for (i in 1:length(pframe.femalemass$timediff))
{
  yvals2[i,] <- vonbertratefunc2(timediff = pframe.femalemass$timediff[i], 
                                  GroupSize.s = pframe.femalemass$GroupSize.s[i],
                                  A = pars.picked2[,1], k = pars.picked2[,2], 
                                  AGS = pars.picked2[,3], kGS = pars.picked2[,4])
}
ci.females <- get_CI(yvals2)
pframe.femalemass <- cbind(pframe.femalemass, ci.females)  

 # male age-related growth
par(mfrow = c(1,2))
plot(weight2 ~ timediff, data = pframe.malemass, type = "n", bty = "l", las = 1, 
     xlab = "Time since first Capture (days)", 
     ylab = "Body Mass (g)", ylim = c(0, 145))
axis(1, at = seq(0, 750, 50), labels = NA)
axis(2, at = seq(10, 150, 10), labels = NA)

with(subset(pframe.malemass, GroupSize.s < -0.5), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[1], alpha.f = 0.2), border=NA))
with(subset(pframe.malemass, GroupSize.s < 0 & GroupSize.s > -0.5), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[2], alpha.f = 0.2), border=NA))
with(subset(pframe.malemass, GroupSize.s > 0), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[3], alpha.f = 0.2), border=NA))
points(weight2 ~ timediff, data = subset(pframe.malemass, GroupSize.s < -0.5), type = "l", lwd = 2, lty =3, col = colpal[1])
points(weight2 ~ timediff, data = subset(pframe.malemass, GroupSize.s< 0 & GroupSize.s > -0.5), type = "l", lwd = 2, lty =2, col = colpal[2])
points(weight2 ~ timediff, data = subset(pframe.malemass, GroupSize.s > 0), type = "l", lwd = 2, lty =1, col = colpal[3])

# female age-related growth
plot(weight2 ~ timediff, data = pframe.femalemass, type = "n", bty = "l", las = 1, 
     xlab = "Time since first Capture (days)", 
     ylab = "Body Mass (g)", ylim = c(0, 145))
axis(1, at = seq(0, 750, 50), labels = NA)
axis(2, at = seq(10, 150, 10), labels = NA)
with(subset(pframe.femalemass, GroupSize.s  < -0.5), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[1], alpha.f = 0.2), border=NA))
with(subset(pframe.femalemass, GroupSize.s > 0 & GroupSize.s < 0.5), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[2], alpha.f = 0.2), border=NA))
with(subset(pframe.femalemass, GroupSize.s > 0.5), polygon(c(timediff, rev(timediff)), c(lwr, rev(upr)), col = adjustcolor(colpal[3], alpha.f = 0.2), border=NA))
points(weight2 ~ timediff, data = subset(pframe.femalemass,  GroupSize.s  < -0.5), type = "l", lwd = 2, lty =3, col = colpal[1])
points(weight2 ~ timediff, data = subset(pframe.femalemass,  GroupSize.s > 0 & GroupSize.s < 0.5), type = "l", lwd = 2, lty =2, col = colpal[2])
points(weight2 ~ timediff, data = subset(pframe.femalemass,  GroupSize.s > 0.5), type = "l", lwd = 2, lty =1, col = colpal[3])
legend(0, 145, legend = c("small", "medium", "large"), lty = c(3,2,1), lwd = 2, col = c(colpal[1], colpal[2], colpal[3]), bty = "n")


# TEETH WIDTH 
#---------------

# first run a basic nls to get an idea of starting parameters for the teeth width models
# female teeth width - baseline

summary(nls(TeethWidth2 ~ A - (A - TeethWidth)*exp(-k*timediff),
            data= maleteeth,
            start=c(A=7,k=0.005))) # A = 6.7586405, k = 0.0033687

vonbert.male.teeth <- nlme(TeethWidth2 ~ A - (A - TeethWidth)*exp(-k*timediff),
                           fixed=A+k~1, random= list(AnimalID = pdDiag(A + k ~ 1)),
                           data= maleteeth,
                           start=c(A=6.75,k=0.005), 
                           na.action = na.omit)
summary(vonbert.male.teeth)   # A = 6.481952; k = 0.004272 

# include group size in the models
# First standardise group size
mean(femaleteeth$GroupSize) ; sd(femaleteeth$GroupSize) # 10.81 # 6.01
mean(maleteeth$GroupSize) ; sd(maleteeth$GroupSize) # 12.52 # 6.01

femaleteeth$GroupSize.s <- as.numeric(scale(femaleteeth$GroupSize))
maleteeth$GroupSize.s <- as.numeric(scale(maleteeth$GroupSize))

# Female model 2
vonbert.female2.teeth <- nlme(TeethWidth2 ~ (A + AGS*GroupSize.s) - ((A + AGS*GroupSize.s) - TeethWidth)*exp(-(k + kGS*GroupSize.s)*timediff),
                              fixed=A+k + AGS + kGS ~1, 
                              random= list(AnimalID = pdDiag(A + k ~ 1)),
                              data= femaleteeth,
                              start=c(A=6,k=0.01, AGS = -0.5, kGS = 0), 
                              na.action= na.omit, 
                              control = nlmeControl(pnlsTol = 0.01, tolerance = 1e-3))
summary(vonbert.female2.teeth)   # A = 5.753846, k = 0.005339, AGS = -0.112547, kGS = 0.000307 


vonbert.male2.teeth <- nlme(TeethWidth2 ~ (A + AGS*GroupSize.s) - ((A + AGS*GroupSize.s) - TeethWidth)*exp(-(k + kGS*GroupSize.s)*timediff),
                            fixed=A+k + AGS + kGS ~1, 
                            random= list(AnimalID = pdDiag(A  + k~ 1)),
                            data= maleteeth,
                            start=c(A=6.75,k=0.005, AGS = -0.5, kGS = 0), 
                            na.action= na.omit)
summary(vonbert.male2.teeth)   # A = 6.478592; k = 0.004281, AGS = -0.047176, kGS = 0.000065

# no real effect on male skeletal growth

AIC(vonbert.female.teeth, vonbert.female2.teeth) %>% 
  mutate(AIC - min(AIC))
AIC(vonbert.male.teeth, vonbert.male2.teeth) %>% 
  mutate(AIC - min(AIC))

summary(vonbert.female2.teeth)
summary(vonbert.male2.teeth)

# PLOT THE TEETH WIDTH RATE OF CHANGE 
#-----------------------------------------------------------

# Plot males and females separately changing the colours according to the predicted group size effect 
#----------------------------------------------

summary(vonbert.male2.teeth)

  plot(rateofincisorgrowth ~ TeethWidth, data = maleteeth, pch = 1, col = adjustcolor("grey20", alpha.f = 0.6), xlab = "Upper incisor Width (g)", ylab = "Rate of change in upper incisor width (mm/day)", bty = "l", las = 1, xlim= c(1.9, 7.5),    ylim = c(-0.0020, 0.015))
  axis(1, at = seq(0, 7.5, 0.5), labels = NA) ; axis(2, at = seq(-0.005, 0.015, 0.0025), labels = NA)
  abline(h = 0, lty =2, col = adjustcolor("grey", alpha.f = 0.6))
  abline(v = 6.472217, lty =1, col = adjustcolor("grey", alpha.f = 0.9))
  
  vonbertratefunc2 <- function(x, A, k, AGs, kGs, GS) {
    time2 = (A + AGs*GS) - ((A + AGs*GS) - x)*exp(-(k + kGs*GS)*180)
    return(time2)
  }
  
  # at 4, 12, and 20 individuals
  small <- (4 - mean(maleteeth$GroupSize))/sd(maleteeth$GroupSize)
  malesmall <- data.frame(time1 = seq(2, 6.51, 0.01),
                          time2 = vonbertratefunc2(x = seq(2, 6.51, 0.01), 
                                                   A = 6.478592, 
                                                   AGs = -0.004281, 
                                                   k = 0.004292, 
                                                   kGs = 0.000065, 
                                                   GS = small)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malesmall, type = "l", lwd = 3, lty = 3, col = colpal[1])
  
  medium <- (12 - mean(maleteeth$GroupSize))/sd(maleteeth$GroupSize)
  malemedium <-  data.frame(time1 = seq(2, 6.47, 0.01),
                            time2 = vonbertratefunc2(x = seq(2, 6.47, 0.01), 
                                                     A = 6.478592, 
                                                     AGs = -0.004281, 
                                                     k = 0.004292, 
                                                     kGs = 0.000065, 
                                                     GS = medium)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malemedium, type = "l", lwd = 3, lty = 2, col = colpal[2])
  
  
  large <- (20 - mean(maleteeth$GroupSize))/sd(maleteeth$GroupSize)
  malelarge <-  data.frame(time1 = seq(2, 6.42, 0.01),
                           time2 = vonbertratefunc2(x = seq(2, 6.42, 0.01), 
                                                    A = 6.478592, 
                                                    AGs = -0.004281, 
                                                    k = 0.004292, 
                                                    kGs = 0.000065, 
                                                    GS = large)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= malelarge, type = "l", lwd = 3, lty = 1, col = colpal[3])
  text(2.5, -0.0015, "Males", cex = 1.2)


# females
#----------

summary(vonbert.female2.teeth)

  plot(rateofincisorgrowth ~ TeethWidth, data = femaleteeth, pch = 4, col = adjustcolor("grey20", alpha.f = 0.6), xlab = "Upper incisor Width (g)", ylab = "Rate of change in upper incisor width (mm/day)", bty = "l", las = 1, xlim= c(1.9, 7.5), ylim = c(-0.0020, 0.015))
  axis(1, at = seq(0, 7.5, 0.5), labels = NA) ; axis(2, at = seq(-0.005, 0.015, 0.0025), labels = NA)
  abline(h = 0, lty =2, col = adjustcolor("grey", alpha.f = 0.6))
  abline(v = 5.754744, lty =1, col = adjustcolor("grey", alpha.f = 0.9))
  
  # at 4, 12, and 20 individuals
  small <- (4 - mean(femaleteeth$GroupSize))/sd(femaleteeth$GroupSize)
  femalesmall <- data.frame(time1 = seq(2, 5.88, 0.01),
                            time2 = vonbertratefunc2(x = seq(2, 5.88, 0.01), 
                                                     A = 5.753846, 
                                                     AGs = -0.112547, 
                                                     k = 0.005339, 
                                                     kGs = 0.000307, 
                                                     GS = small)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalesmall, type = "l", lwd = 3, lty = 3, col = colpal[1])
  
  medium <- (12 - mean(femaleteeth$GroupSize))/sd(femaleteeth$GroupSize)
  femalemedium <-  data.frame(time1 = seq(2, 5.75, 0.01),
                              time2 = vonbertratefunc2(x = seq(2, 5.75, 0.01), 
                                                       A = 5.754744, 
                                                       AGs = -0.112128, 
                                                       k = 0.005336, 
                                                       kGs = 0.000305, 
                                                       GS = medium)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalemedium, type = "l", lwd = 3, lty = 2, col = colpal[2])
  
  
  large <- (20 - mean(femaleteeth$GroupSize))/sd(femaleteeth$GroupSize)
  femalelarge <-  data.frame(time1 = seq(2, 5.62, 0.01),
                             time2 = vonbertratefunc2(x = seq(2, 5.62, 0.01), 
                                                      A = 5.754744, 
                                                      AGs = -0.112128, 
                                                      k = 0.005336, 
                                                      kGs = 0.000305, 
                                                      GS = large)) %>% 
    mutate(rate = (time2 - time1)/180)
  points(rate ~ time1, data= femalelarge, type = "l", lwd = 3, lty = 1, col = colpal[3])
  text(2.5, -0.0015, "Females", cex = 1.2)


#------------------------------------------------------------------------------------------  
  

# Body mass graphs indicate that adults in larger groups are of lower asymptotic mass. 
# Just to check this trend in a more conventional manner we will visuliase the adult mass of all individuals 
  # that were first caught at less than 80g or 100g (~ 1 year of age in either sex)
  
malesunder1yr <- malebodymass %>% 
    mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y")), 
           DateofRecapture = as.Date(strptime(DateofRecapture, format = "%d/%m/%Y"))) %>% 
    arrange(AnimalID, DateofCapture) %>% 
    group_by(AnimalID) %>% 
    slice(1) %>% 
    ungroup() %>% 
    filter(weight < 101) %>% 
    dplyr::select(AnimalID, DateofCapture) %>% 
    rename(FirstCapture = DateofCapture)

# now get the weights of all these individuals at least 1 year beyond that first weight (guarenteeing they are an adult by then) 
# do in two stages
  
massover1yr <- malebodymass %>% 
  dplyr::select(-DateofRecapture, -weight2, -GroupSizeRecapture, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
  right_join(malesunder1yr, by = "AnimalID") %>%
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
  mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
  filter(timesincefirstcapture > 365) %>% 
  bind_rows(malebodymass %>% 
  dplyr::select(-DateofCapture, -weight, -GroupSize, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
  rename(DateofCapture = DateofRecapture, weight = weight2, GroupSize = GroupSizeRecapture) %>% 
  right_join(malesunder1yr, by = "AnimalID") %>%
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
  mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
  filter(timesincefirstcapture > 365)) %>% 
  distinct()

# males
quantile(massover1yr$GroupSize, probs = c(0.333, 0.666))
massover1yr <- massover1yr %>% 
  mutate(GroupSizeCat = as.factor(case_when(GroupSize %in% 1:9 ~ "Small", 
                                  GroupSize %in% 10:16 ~ "Medium", 
                                  GroupSize > 16 ~ "Large")))
massover1yr$GroupSizeCat <- factor(massover1yr$GroupSizeCat, levels = c("Small", "Medium", "Large"))

par(mfrow = c(1,2))
plot(weight ~ as.factor(GroupSizeCat), data = massover1yr, xlab = "Group Size Category", ylab = "Weight (g)", las = 1, 
     col = c(colpal[1],colpal[2], colpal[3]), bty = "n", varwidth = TRUE, boxwex =0.6, ylim = c(60, 202))
table(massover1yr$GroupSizeCat)
text(c(1, 2, 3), c(202, 202, 202), c(47, 56, 49), cex = 0.9)

#####

femalesunder1yr <- femalebodymass %>% 
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y")), 
         DateofRecapture = as.Date(strptime(DateofRecapture, format = "%d/%m/%Y"))) %>% 
  arrange(AnimalID, DateofCapture) %>% 
  group_by(AnimalID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(weight < 81) %>% 
  dplyr::select(AnimalID, DateofCapture) %>% 
  rename(FirstCapture = DateofCapture)


# now get the weights of all these individuals at least 1 year beyond that first weight (guarenteeing they are an adult by then) 
# do in two stages

massover1yrfem <-   femalebodymass %>% 
  dplyr::select(-DateofRecapture, -weight2, -GroupSizeRecapture, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
  right_join(femalesunder1yr, by = "AnimalID") %>%
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
  mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
  filter(timesincefirstcapture > 365) %>% 
  bind_rows(femalebodymass %>% 
              dplyr::select(-DateofCapture, -weight, -GroupSize, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
              rename(DateofCapture = DateofRecapture, weight = weight2, GroupSize = GroupSizeRecapture) %>% 
              right_join(femalesunder1yr, by = "AnimalID") %>%
              mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
              mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
              filter(timesincefirstcapture > 365)) %>% 
  distinct()

# females
quantile(massover1yrfem$GroupSize, probs = c(0.333, 0.666)) # 8. 15
massover1yrfem <- massover1yrfem %>% 
  mutate(GroupSizeCat = as.factor(case_when(GroupSize %in% 1:7 ~ "Small", 
                                            GroupSize %in% 8:15 ~ "Medium", 
                                            GroupSize > 15 ~ "Large")))
massover1yrfem$GroupSizeCat <- factor(massover1yrfem$GroupSizeCat, levels = c("Small", "Medium", "Large"))


plot(weight ~ as.factor(GroupSizeCat), data = massover1yrfem, xlab = "Group Size Category", ylab = "Weight (g)", las = 1, 
     col = c(colpal[1],colpal[2], colpal[3]), bty = "n", varwidth = TRUE, boxwex =0.6, ylim = c(60, 170))
table(massover1yrfem$GroupSizeCat)
text(c(1, 2, 3), c(172, 172, 172), c(37, 38, 36), cex = 0.9)


# Include a term for the early life group size of each individual as well and plot. 

# Females: Group size in the first year of their life
groupsizeunder1year <- femalebodymass %>% 
  dplyr::select(-DateofRecapture, -weight2, -GroupSizeRecapture, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
  right_join(femalesunder1yr, by = "AnimalID") %>%
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
  mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
  filter(weight < 81) %>% 
  group_by(AnimalID) %>% 
  summarise(avgGroupSize.1yr = mean(GroupSize),   
            FirstGroupSize = GroupSize[which(DateofCapture == min(DateofCapture))], 
            n.GroupSizeMeasures = n())


# Males: 
groupsizeunder1year.males <- malebodymass %>% 
  dplyr::select(-DateofRecapture, -weight2, -GroupSizeRecapture, -QueenID, -timediff, -rateofbodymassgrowth, -GroupSize.s) %>% 
  right_join(malesunder1yr, by = "AnimalID") %>%
  mutate(DateofCapture = as.Date(strptime(DateofCapture, format = "%d/%m/%Y"))) %>% 
  mutate(timesincefirstcapture = as.numeric(difftime(DateofCapture, FirstCapture, units = c("days")))) %>% 
  filter(weight < 101) %>% 
  group_by(AnimalID) %>% 
  summarise(avgGroupSize.1yr = mean(GroupSize),   
            FirstGroupSize = GroupSize[which(DateofCapture == min(DateofCapture))], 
            n.GroupSizeMeasures = n())

massover1yrfem <- massover1yrfem %>% 
  left_join(groupsizeunder1year)

massover1yr <- massover1yr %>% 
  left_join(groupsizeunder1year.males)

# Statistically test these effects #
library(lme4) ; library(multcomp) 

# In males: 
#-----------

# include time of year
massover1yr <- massover1yr %>% 
  mutate(month = lubridate::month(DateofCapture), 
         MonthlyQuarter = as.factor(case_when(month %in% 1:3 ~ "first", 
                                    month %in% 4:6 ~ "second", 
                                    month %in% 7:9 ~ "third", 
                                    month %in% 10:12 ~ "fourth")))


# how correlated are average group size below one year and group size in adulthood
plot(GroupSize ~ avgGroupSize.1yr, data = massover1yr) # males- not that correlated 
abline(lm(GroupSize ~ avgGroupSize.1yr, data = massover1yr), col = "red")
cor.test(massover1yr$GroupSize, massover1yr$avgGroupSize.1yr)

massover1yr$avgGroupSize.1yr.s <- as.numeric(scale(massover1yr$avgGroupSize.1yr))
mod1 <- lme4::lmer(weight ~ GroupSizeCat + MonthlyQuarter + avgGroupSize.1yr.s + (1|GroupID/AnimalID), data = massover1yr) # very small differences in reality in both sexes in adulthood

summary(mod1)
anova(mod1)
plot(mod1)

anova(mod1, update(mod1, ~. -GroupSizeCat))
anova(mod1, update(mod1, ~. -MonthlyQuarter))
anova(mod1, update(mod1, ~. -avgGroupSize.1yr.s))

summary(glht(mod1, linfct = mcp(GroupSizeCat = "Tukey")), test = adjusted("holm"))

# plot the marginal effects
p1 <- ggeffects::ggpredict(mod1, terms = "GroupSizeCat") %>% 
  data.frame() %>% 
  rename(GroupSizeCat = x, weight = predicted)

plot1 <- ggplot(p1, aes(GroupSizeCat, weight)) + 
  geom_errorbar(aes(ymin=weight-1.96*std.error, ymax=weight+1.96*std.error), width=.2, size = 1.2) + 
  geom_jitter(data = massover1yr, aes(col = GroupSizeCat),  width = 0.1, size = 2, alpha = 0.4) + 
  scale_colour_manual(values = c("#56B4E9","#009E73", "#0072B2")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") + 
  ylab("Body Mass, g") + 
  xlab("Group Size Category") 

summary(mod1)

p2 <- ggeffects::ggpredict(mod1, terms = "avgGroupSize.1yr.s")

plot2 <- plot(p2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") + 
  ylab("Body Mass, g") + 
  xlab("Average group size in first year of life") 


# In Females: Female mass in adulthood
#--------------------------------------

massover1yrfem <- massover1yrfem %>% 
  mutate(month = lubridate::month(DateofCapture), 
         MonthlyQuarter = as.factor(case_when(month %in% 1:3 ~ "first", 
                                    month %in% 4:6 ~ "second", 
                                    month %in% 7:9 ~ "third", 
                                    month %in% 10:12 ~ "fourth")))


massover1yrfem$avgGroupSize.1yr.s <- as.numeric(scale(massover1yrfem$avgGroupSize.1yr))
mod2 <- lme4::lmer(weight ~ GroupSizeCat + MonthlyQuarter  + avgGroupSize.1yr.s + (1|GroupID/AnimalID), data = massover1yrfem) # very small differences in reality in both sexes in adulthood
summary(mod2)
anova(mod2)
plot(mod2)

anova(mod2, update(mod2, ~. -GroupSizeCat))
anova(mod2, update(mod2, ~. -MonthlyQuarter))
anova(mod2, update(mod2, ~. -avgGroupSize.1yr.s))

summary(glht(mod2, linfct = multcomp::mcp(GroupSizeCat = "Tukey")), test = adjusted("holm"))

# plot the marginal effects
p3 <- ggeffects::ggpredict(mod2, terms = "GroupSizeCat") %>% 
  data.frame() %>% 
  rename(GroupSizeCat = x, weight = predicted)

plot3 <- ggplot(p3, aes(GroupSizeCat, weight)) + 
  geom_errorbar(aes(ymin=weight-1.96*std.error, ymax=weight+1.96*std.error), width=.2, size = 1.2) + 
  geom_jitter(data = massover1yrfem, aes(col = GroupSizeCat),  width = 0.1, size = 2, alpha = 0.4) + 
  scale_colour_manual(values = c("#56B4E9","#009E73", "#0072B2")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") + 
  ylab("Body Mass, g") + 
  xlab("Group Size Category") 


p3 <- ggeffects::ggpredict(mod2, terms = "avgGroupSize.1yr.s")

plot4 <- plot(p3) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.position = "none") + 
  ylab("Body Mass, g") + 
  xlab("Average group size in first year of life") 
plot4

#####-----------------  END   ------------------####
