#----------------------------------------------------------------
# No obvious benefit of group size in wild Damaraland mole-rats Fukomys damarensis
#
# i- "Status-related survivoship in females" 
#
# R script
# Authors: Jack Thorley, Hanna Bensch, Markus Zottl
# Contact: jackthorley1@gmail.com
#----------------------------------------------------------------

library(tidyverse) ;  library(msm) ; library(msmtools)

# In this script we explore the life history trajectories of females across all captures.  
# We  focus on all females (irrespective of when they were first captured), and build a single multi-state model in order to investigate differences in 'survivorship' the between three states: non-breeders in their group, single females, and breeding females. 
# Technically, because we cannot know that individuals who disappear are dead, we are really estimating disappearance. 

# read in the appropriate data set (will need to set your own working directory here)
females <- read.csv("Fates_AllFemales.csv", header = TRUE) # Females(first captured < 80g)

# first need to incorporate the plotting windows 
trappingwindows.krr <- data.frame(TrappingWindow = 1:13, PeriodStart = as.Date(c("2013-11-17", "2014-07-27", "2015-02-10", "2015-09-16", "2016-02-04", "2016-07-25", "2017-01-16", "2017-07-31", "2018-03-17", "2018-09-14", "2019-03-09", "2019-09-09", "2020-03-05")), PeriodEnd = as.Date(c( "2014-05-21","2014-11-24", "2015-04-14", "2015-10-22", "2016-06-06", "2016-09-12", "2017-05-30", "2017-10-26", "2018-05-30", "2018-12-04", "2019-04-24", "2019-11-30", "2020-05-27")))

trappingwindows.lonely <- data.frame(TrappingWindow = 1:7, PeriodStart = as.Date(c("2013-09-23", "2014-03-01", "2014-09-16", "2015-01-29", "2016-02-08", "2016-09-09", "2017-05-30")), PeriodEnd = as.Date(c("2013-10-06", "2014-07-03", "2014-10-27", "2015-04-02", "2016-05-12", "2016-10-13", "2017-08-02")))

# The plot for young females from "KRR"
females$CaptureStart <- as.Date(strptime(females$CaptureStart, format = "%d/%m/%Y"))
females$firstcapturedate <- as.Date(strptime(females$firstcapturedate, format = "%d/%m/%Y"))

r1 <- females %>% 
  group_by(AnimalID) %>% 
  mutate(firstcapturedate = min(CaptureStart)) %>% 
  ungroup() %>% 
  dplyr::select(AnimalID, firstcapturedate) %>% 
  distinct() %>% 
  arrange(firstcapturedate)

r2 <- rle(as.character(r1$AnimalID))
r1$id <- rep(seq_along(r2$lengths), r2$lengths)
r1 <- select(r1, -firstcapturedate)
  
females <- females %>% 
  left_join(r1) %>% 
  arrange(id) %>% 
  mutate(id = as.factor(id), statecode = as.factor(statecode))

females$statecode <- factor(females$statecode, levels = c("femalenonbreeder", "single", "inheritedbreeder", "breeder",  "Dead"))

female.fates <- ggplot(females, aes(x = CaptureStart, y = id, group = id, fill = statecode, shape = statecode))  +
  geom_vline(xintercept = seq(365/2, 2200, 365/2), colour="grey", linetype = "longdash") +
  scale_x_date(date_breaks = "6 month", date_labels = "%b %Y") +
  geom_line(alpha = 0.7, col = "black") + 
  xlab('Capture Date') + 
  ylab('Animal ID') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6), 
        axis.text.y = element_blank()) + 
  geom_rect(data=trappingwindows.krr, inherit.aes=FALSE,
            aes(xmin=PeriodStart,xmax=PeriodEnd,ymin=0,ymax=364), alpha = 0.2, fill = "lightblue")  + 
  geom_rect(data=trappingwindows.lonely, inherit.aes=FALSE,
            aes(xmin=PeriodStart,xmax=PeriodEnd,ymin=0,ymax=364), alpha = 0.2, fill = "lightgreen")  + 
  geom_point(size  = 4, alpha = 0.6) + 
  scale_shape_manual(values = c(21, 21, 21, 21, 21), labels = c("Female nonbreeder", "Single female (dispersed)", "Inherited breeding position", "Breeder", "Dead/Disappeared")) +
  scale_fill_manual(values  = c("gold", "darkorange", "dodgerblue", "forestgreen", "red"), guide = FALSE) + 
  scale_colour_manual(values = c("black", "black", "black", "black")) + 
  guides(shape = guide_legend(override.aes = list(fill = c("gold", "darkorange", "dodgerblue", "forestgreen", "red")))) +
  ggtitle("Female mole-rat captures")

annotation <- females %>% 
  filter(censor == 99) %>% 
  mutate(label = "X")

female.fates <- female.fates +  
  geom_text(data=annotation, aes(CaptureStart, y = id, label=label), 
            color="black", 
            size=5, fontface="bold")
female.fates

# Note that two uncommon transitions occur here: 
# L40F002 and KR13F017 both became breeders but then subsequently became single (the group collapsed around them)


#----------------------------- Fitting a global multistate model

females <- arrange(females, AnimalID, CaptureStart)

statetable.msm(statecode, AnimalID, data = females)
unique(filter(females, statecode == "inheritedbreeder")$AnimalID)

# generate a state variable
# because very few individuals in their natal group inherit the breeding position, inherited breeders and out-of-group breeders are collapsed into a single state (state '3')
# then we can allow for female non-breeders to become a breeder through two routes (either directly via 1, or indirectly via 2)

females <- females %>% 
  mutate(state = case_when(statecode == "femalenonbreeder" ~ 1, 
                           statecode == "single" ~ 2,
                           statecode == "inheritedbreeder" ~ 3,
                           statecode == "breeder" ~ 3, 
                           statecode == "Dead" ~ 4))
females$state[which(!is.na(females$censor))] <- 99

# In a souple of places they go from breeder to single individual, and breeder to female non-breeder
statetable.msm(state, AnimalID, data = females)

#  set up the q matrix for the transitions 
# 4-by-4 matrix 
# all diagonals zero, 
females.q <- rbind(c(0, 0.166, 0.166, 0.166),  # from non-breeder to all fates
                   c(0, 0, 0.25, 0.25),            # from single to breeder or dead
                   c(0, 0, 0, 0.5),              # from breeder death
                   c(0, 0, 0 , 0))
rownames(females.q) <- colnames(females.q) <- c("femalenonbreeder", "single", "breeder","dead")

# simple model (assuming either single or paired states can be censored)
# and don't specify death equals as the timing of the final state is also not known with certainty
females.msm <- msm(state ~ Time.years, subject = AnimalID, censor = 99, censor.states = c(1,2,3,4),  
                                data = females, qmatrix = females.q, 
                                na.action = na.fail)
pmatrix.msm(females.msm, t = 1, ci = "normal") # Transition intensities over a period of 1 

# Probability that individual that enters into each state is dead one year later 
survivalprobabilities <- pmatrix.msm(females.msm, t = 2, ci = "normal")[1:3,4]

#survivalprobabilities.bootstrap <- pmatrix.msm(females.msm, t = 1, ci = "bootstrap")[1:3,4]

# Compare the ratio of transition intensities for death at each state: 
# Specifically, the ratio of two entries of the transition intensity matrix at a given set of covariate values, together with a confidence interval, estimated assuming normality on the log scale and using the delta method.
qratio.msm(females.msm, ind1 = c(2,4), ind2 = c(3,4)) # single female versus breeder, not sig
qratio.msm(females.msm, ind1 = c(1,4), ind2 = c(2,4)) # non-breeder versus single female, significant
qratio.msm(females.msm, ind1 = c(1,4), ind2 = c(3,4)) # non-breeder versus breeding female, significant

# Run some simple plots of summary output. 
par(mfrow = c(1,1))
plot.msm(females.msm, las = 1, xlab = "Time, years", ylab = "Survival Probability", cex.lab = 1.2, cex.axis = 1.3)  # note here that survival means not entering absorbing state, so this nicely demonstrates already


par(mfrow = c(2,2))
plot.survfit.msm(females.msm, from = 1, las = 1, main = "In-group non-breeders", interp = "midpoint", lwd = 3, col = "gold", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1)
plot.survfit.msm(females.msm, from = 2, las = 1, main = "Single females", interp = "midpoint", lwd = 3, col = "darkorange", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1)  
plot.survfit.msm(females.msm, from = 3, las = 1, main = "Breeding females", interp = "midpoint", lwd = 3, col = "forestgreen", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1) 

# incorporate confidence intervals in the fitted survivorship
  par(mfrow = c(2,2))
  plot.survfit.msm(females.msm, from = 1, las = 1, main = "In-group non-breeders", interp = "midpoint", lwd = 3, col = "gold", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1, mark.time= FALSE, ci = "normal", col.ci = "gold")
  plot.survfit.msm(females.msm, from = 2, las = 1, main = "Single females", interp = "midpoint", lwd = 3, col = "darkorange", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1, mark.time= FALSE, ci = "normal", col.ci = "darkorange") 
  plot.survfit.msm(females.msm, from = 3, las = 1, main = "Breeding females", interp = "midpoint", lwd = 3, col = "forestgreen", col.surv = adjustcolor("blue", alpha.f = 0.5), cex.lab = 1.1, cex.axis = 1.1, mark.time= FALSE, ci = "normal", col.ci = "forestgreen")


# Lastly plot the annual survival estimate
survivalprobabilities <- data.frame(pmatrix.msm(females.msm, t = 1, ci = "normal")[1:3,4])
df <- data.frame(estimate = survivalprobabilities[1], lower = survivalprobabilities[2], upper = survivalprobabilities[3], 
                 class = c("nonbreeder", "single", "breeder"), num = 1:3)

par(mfrow=c(1,1))
plot(estimate ~ num, data = df, las = 1, type  = "n", xlim = c(0, 4), ylim= c(0, 0.5), xaxt = "n", xlab = NA, ylab = "Annual Probability of Disappearing", bty = "l")
for(i in 1:3){
  lines(c(i, i), c(df$lower[i], df$upper[i]), lwd = 2, col = c("gold", "darkorange", "forestgreen")[i])
}
points(estimate ~ num, data = df, las = 1, type  = "p", col = c("gold", "darkorange", "forestgreen"), cex = 2, pch = 16)
points(estimate ~ num, data = df, las = 1, type  = "p", cex = 2, pch = 1)
axis(1, at = 1:3, labels = NA)

#####-----------------  END   ------------------####