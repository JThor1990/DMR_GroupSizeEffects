#----------------------------------------------------------------
# Small group size does not impair reproduction and survival in social Damaraland mole-rats
#
# v- "Analysing the timing of natal dispersal (natal philopatry) for both sexes 
#
# R script
# Authors: Jack Thorley, Hanna Bensch, Markus Zöttl
# Contact: jack.thorley1@gmail.com
#----------------------------------------------------------------

library(tidyverse) ;  library(msm) ; library(msmtools)


# call in the data
femalephilopatry.df <- read.csv("FieldMR_Femalephilopatry.csv", header = TRUE)
malephilopatry.df <- read.csv("FieldMR_Malephilopatry.csv", header = TRUE)

philopatry.df <- rbind(femalephilopatry.df, malephilopatry.df)


# Model the sexes together 
#---------------------------------------------------------

# Make sure they have the same column names so that they can be bound together
#names(malephilopatry.df)[which(!(names(malephilopatry.df) %in% names(femalephilopatry.df)))]
#names(femalephilopatry.df)[which(!(names(femalephilopatry.df) %in% names(malephilopatry.df)))]

#malephilopatry.df <- malephilopatry.df[, -which(!(names(malephilopatry.df) %in% names(femalephilopatry.df)))]
#femalephilopatry.df <- femalephilopatry.df[, -which(!(names(femalephilopatry.df) %in% names(malephilopatry.df)))]

#malephilopatry.df <- malephilopatry.df[,names(femalephilopatry.df)]

statetable.msm(state, AnimalID, data = philopatry.df)

# set up the q matrix for the transitions (assume that individuals can transition from non-breeder to gone or known out of group)
philopatry.q <- rbind(c(0, 0.5), 
                      c(0, 0)) # assuming about 1 year in state 1 before disappearing

philopatry.df$Sex <- as.factor(philopatry.df$Sex)

rownames(philopatry.q) <- colnames(philopatry.q) <- c("ingroup", "Disappeared/Dispersed")

philopatry.msm <- msm(state ~ Time.years, subject = AnimalID, censor = 99, censor.states = c(1),
                      covariates = ~ Sex, data = philopatry.df, qmatrix = philopatry.q, na.action = na.fail, 
                      control = list(fnscale = 4000, maxit = 10000))
philopatry.msm
summary(philopatry.msm)
hazard.msm(philopatry.msm) # n.s.

#-------------------------------------------------

# Now, refit the model with the two sexes with a weight term and a rainfall term included to see whether: 
# i) heavier individuals are more likely to disappear. 
# ii) rainfall predicts disappearance

# get rid of individuals where weight is an NA that aren't a final event
weightless <- philopatry.df %>% 
  filter(is.na(weight), statecode != "Disappeared", is.na(censor)) %>% 
  .$AnimalID

philopatry.df2 <- philopatry.df %>% 
  filter(!(AnimalID %in% weightless)) 

# easiest to remove the NAs before rebinding
philopatry.df2 <- philopatry.df2 %>% 
  filter(!is.na(weight)) %>% 
  group_by(Sex) %>% 
  mutate(ScaledWeight = as.numeric(scale(weight))) %>% 
  bind_rows(philopatry.df2 %>% 
              filter(is.na(weight))) %>% 
  arrange(AnimalID, CaptureStart) %>% 
  data.frame()

philopatry.q2 <- rbind(c(0, 0.5), 
                       c(0, 0)) # set up the q matrix
rownames(philopatry.q2) <- colnames(philopatry.q) <- c("ingroup", "Disappeared/Dispersed")

philopatry.df2$Rainfall.s <- as.numeric(scale(philopatry.df2$Rainfall))

# model
philopatry.msm.2 <- msm(state ~ Time.years, subject = AnimalID, censor = 99, censor.states = c(1),
                        covariates = ~ Sex + ScaledWeight + Rainfall.s, data = philopatry.df2, 
                        qmatrix = philopatry.q, na.action = na.fail, 
                        control = list(fnscale = 4000, maxit = 10000))
philopatry.msm.2
summary(philopatry.msm.2)
hazard.msm(philopatry.msm.2)


# plot the outputs
par(mfrow = c(1,1))
plot.survfit.msm(philopatry.msm.2, from = 1, las = 1) # overall

par(mfrow = c(1,2))
plot.survfit.msm(philopatry.msm.2, from = 1, las = 1, 
                 covariates = list(Sex = "Female", "Rainfall.s" = 0, "ScaledWeight"  = 0), 
                 interp = "midpoint", xlab = "Time since first captured, years", 
                 ylab = "Probability of presence in natal group", range = c(0,4)) 
plot.survfit.msm(philopatry.msm, from = 1, las = 1, 
                 covariates = list(Sex = "Male", "Rainfall.s" = 0, "ScaledWeight"  = 0),
                 interp = "midpoint", xlab = "Time since first captured, years", 
                 ylab = "Probability of presence in natal group", range = c(0,4))  

# note that fits a different curve but maintains the empirical data so need to extract this separately for males and females from above and combine the plots together

# get the prevalences to compare Females to Males
plot.prevalence.msm(philopatry.msm.2, xlab =  "Time since first captured, years")

philopatry.prevalence.msm.females <- msm:::prevalence.msm(philopatry.msm.2, ci = "normal", 
                                                          covariates = list(Sex = "Female", 
                                                                            Rainfall.s = 0, 
                                                                            ScaledWeight = 0))
philopatry.prevalence.msm.males <- msm:::prevalence.msm(philopatry.msm.2, ci = "normal", 
                                                        covariates = list(Sex = "Male", 
                                                                          Rainfall.s = 0, 
                                                                          ScaledWeight = 0))

#  tidy into a easier format
# females, observed prevalence
femalephilopatry.observed <- philopatry.prevalence.msm.females[[3]] %>% 
  data.frame() %>% 
  mutate(Time = as.numeric(row.names(.))) %>% 
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "observedprevalence")

# females, expected percentage/prevalence - mean
femalephilopatry.expected <- philopatry.prevalence.msm.females$`Expected percentages`$estimates %>% 
  data.frame() %>% 
  mutate(Time = as.numeric(row.names(.))) %>% 
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "expectedprevalence")

#  females, expected percentage/prevalence - l95%CI
femalephilopatry.expected.lci <- philopatry.prevalence.msm.females$`Expected percentages`$ci[,,1]%>% 
  data.frame() %>% 
  mutate(Time = unique(femalephilopatry.expected$Time)) %>%  
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "l95ci") 

#  females, expected percentage/prevalence - u95%CI
femalephilopatry.expected.uci <- philopatry.prevalence.msm.females$`Expected percentages`$ci[,,2]%>% 
  data.frame() %>% 
  mutate(Time = unique(femalephilopatry.expected$Time)) %>%  
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "u95ci") 

femalephilopatry.expected <- femalephilopatry.expected %>% # put together
  left_join(femalephilopatry.expected.lci) %>% 
  left_join(femalephilopatry.expected.uci) %>% 
  filter(state== "ingroup")

# males, observed prevalence
malephilopatry.observed <- philopatry.prevalence.msm.males[[3]] %>% 
  data.frame() %>% 
  mutate(Time = as.numeric(row.names(.))) %>% 
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "observedprevalence")

# males, expected percentage/prevalence - mean
malephilopatry.expected <- philopatry.prevalence.msm.males$`Expected percentages`$estimates %>% 
  data.frame() %>% 
  mutate(Time = as.numeric(row.names(.))) %>% 
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "expectedprevalence")

# males, expected percentage/prevalence - l95%CI
malephilopatry.expected.lci <- philopatry.prevalence.msm.males$`Expected percentages`$ci[,,1]%>% 
  data.frame() %>% 
  mutate(Time = unique(malephilopatry.expected$Time)) %>%  
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "l95ci") 

# males, expected percentage/prevalence - U95%CI
malephilopatry.expected.uci <- philopatry.prevalence.msm.males$`Expected percentages`$ci[,,2]%>% 
  data.frame() %>% 
  mutate(Time = unique(malephilopatry.expected$Time)) %>%  
  rename(ingroup = names(.)[1], outgroup = names(.)[2]) %>% 
  pivot_longer(-Time, names_to = "state", values_to = "u95ci") 

malephilopatry.expected <- malephilopatry.expected %>%  # put together 
  left_join(malephilopatry.expected.lci) %>% 
  left_join(malephilopatry.expected.uci) %>% 
  filter(state== "ingroup") %>% 
  data.frame()

# --> plots, predicted philopatry of males and females with CI 
plot(expectedprevalence ~ Time, data = subset(femalephilopatry.expected, state == "ingroup"), type = "n", las = 1, bty = "l", ylab= "Probability of presence in natal group (%)", xlab = 'Time since first captured \n as non-breeding male, years', ylim = c(0, 101), xlim = c(0, 5))
with(femalephilopatry.expected, 
     polygon(c(Time, rev(Time)), c(l95ci, rev(u95ci)), 
             col = adjustcolor("#1F9E89FF", alpha = 0.2), border = NA))
with(malephilopatry.expected, 
     polygon(c(Time, rev(Time)), c(l95ci, rev(u95ci)), 
             col = adjustcolor("black", alpha = 0.1), border = NA))
points(expectedprevalence ~ Time, 
       data = femalephilopatry.expected, type = "l", col = "#1F9E89FF", lwd = 2) # expected
points(expectedprevalence ~ Time, data = subset(malephilopatry.expected, state == "ingroup"), type = "l", col = "black", lwd = 2) # expected
axis(1, at = seq(0, 5, 0.5), labels = NA)
axis(2, at = seq(0, 100, 10), labels = NA)
legend("topright", legend = c("Female", "Male"), lwd = c(2,2), bty = "n", col = c("#1F9E89FF", "black"))

# mean sojourn times for males and females in their natal group
sojourn.msm(philopatry.msm.2, covariates = list(Sex = "Female", 
                              Rainfall.s = 0, 
                              ScaledWeight = 0))

sojourn.msm(philopatry.msm.2, covariates = list(Sex = "Male", 
                                                Rainfall.s = 0, 
                                                ScaledWeight = 0))
