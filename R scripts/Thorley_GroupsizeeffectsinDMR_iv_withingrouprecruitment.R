#----------------------------------------------------------------
# No obvious benefit of group size in wild Damaraland mole-rats Fukomys damarensis
#
# iv- "Within-group recruitment rate" 
#
# R script
# Authors: Jack Thorley, Hanna Bensch, Markus Zottl
# Contact: jackthorley1@gmail.com
#----------------------------------------------------------------

library(glmmTMB) ; library(tidyverse)

# In this script we explore the within-group recruitment rate in a longitudinal analysis of groups, and following experimental pairings. 

# load in the recruitment datasets 
longitudinal <- read.csv("FieldMR_Recruitment_Longitudinal.csv", header = TRUE)  # "FieldMR_Recruitment_Longitudinal.csv"
experimental <- read.csv("FieldMR_Recruitment_Experimental.csv", header = TRUE)  # "FieldMR_Recruitment_Experimental.csv"

# Model within-group recruitment in the longitudinal data
#------------------------------------------------------------

#  Scale the various continous variables
longitudinal$GroupSize.s <- as.numeric(scale(longitudinal$GroupSizeAtFirstCapture))
longitudinal$QueenWeight.s <- as.numeric(scale(longitudinal$QueenWeight))
longitudinal$Rainfall1.s <- as.numeric(scale(longitudinal$Rainfall1)) # 1) sum rainfall to the year before the start of the capture period
longitudinal$Rainfall2.s <- as.numeric(scale(longitudinal$Rainfall2)) # 2) Geometric mean rainfall in the year before the start 
longitudinal$Rainfall3.s <- as.numeric(scale(longitudinal$Rainfall3)) # 3) Arithmetic mean rainfall in the year before the start 
longitudinal$Rainfall4.s <- as.numeric(scale(longitudinal$Rainfall4)) # 4) sum within the trapping interval
longitudinal$Rainfall5.s <- as.numeric(scale(longitudinal$Rainfall5)) # 5) geometric mean monthly in the trapping interval
longitudinal$Rainfall6.s <- as.numeric(scale(longitudinal$Rainfall6)) # 6) Arithmetic mean rainfall in the year before the start 

# fit a poisson and negative binomial to a baseline model first to see what is going on

mod1 <- glmmTMB(WithinGroupRecruits ~ GroupSize.s + QueenWeight.s + offset(log(TimeToNextCap))
                + (1|GroupID), 
                family = "poisson", 
                data = longitudinal)

mod2 <- update(mod1, ~., family = "nbinom2")

AIC(mod1, mod2) # Here the poisson is a better fit 

# overdispersed?
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(mod1) 

# check residuals of the Poisson models
simulationOutput.mod1 <- DHARMa::simulateResiduals(fittedModel = mod1, n = 250)
plot(simulationOutput.mod1) # look alright
hist(simulationOutput.mod1)


# random effects significant?
anova(mod1, update(mod1, ~. -(1|GroupID))) # groupID significant at alpha = 0.05 in the Poisson

# Okay now to this  model in each case to the different rainfall terms and see which is 'best'
mod1.1 <- update(mod1, ~. + Rainfall1.s)
mod1.2 <- update(mod1, ~. + Rainfall2.s) # best fitting rainfall term
mod1.3 <- update(mod1, ~. + Rainfall3.s)
mod1.4 <- update(mod1, ~. + Rainfall4.s)
mod1.5 <- update(mod1, ~. + Rainfall5.s)
mod1.6 <- update(mod1, ~. + Rainfall6.s)

summary(mod1.2)
AIC(mod1.1, mod1.2, mod1.3, mod1.4, mod1.5, mod1.6)
# total rainfall in the preceding year is apparently 'best' (most correlated) rainfall metric

newdat <- data.frame(GroupSize.s = (seq(2, 26, 0.5) - mean(longitudinal$GroupSizeAtFirstCapture))/sd(longitudinal$GroupSizeAtFirstCapture),   QueenWeight.s = 0, TimeToNextCap = 180, NewRecruits = 0, Rainfall2.s = 9, GroupID = NA)

preds <- predict(mod1,newdat, re.form=NA, se.fit = TRUE) 
newdat$pred <- exp(preds$fit)
newdat$ulimit <- exp(preds$fit + qnorm(0.975)*(preds$se.fit)) # needs to be on the correct scale when getting them
newdat$llimit <- exp(preds$fit - qnorm(0.975)*(preds$se.fit)) 
newdat$GroupSize <- seq(2, 26, 0.5)

plot(withingrouprecruitmentrate.6mo ~ GroupSizeAtFirstCapture, data = longitudinal, las = 1, type = "n",
     ylab = "Recruitment rate, within-group recruits per 6 months", xlab = "Group Size, t", 
     xlim = c(0, 26), 
     bty = "l")
axis(1, at = seq(0, 26, 1), label = NA)
with(newdat, polygon(c(GroupSize, rev(GroupSize)), c(ulimit, rev(llimit)), col = adjustcolor("darkorange", alpha.f = 0.3), border = NA))
with(newdat, points(pred ~ GroupSize, type = "l", col = "darkorange", lwd = 3))
points(withingrouprecruitmentrate.6mo ~ GroupSizeAtFirstCapture, data = longitudinal, col =  adjustcolor("black", alpha.f = 0.9), cex = 1.2)

mean(longitudinal$TimeToNextCap) ; sd(longitudinal$withingrouprecruitmentrate.6mo)

# a very weak group size effect
summary(mod1.2)

# re-run the model with other forms of group size term
summary(update(mod1.2, ~.-GroupSize.s + AvgGroupSize)) 
summary(update(mod1.2, ~.-GroupSize.s + AvgGroupSizeWoutRecruits)) # doesn't make much difference. 


# Experimental investigation of recruitment
#--------------------------------------------

# Compare the two 'tratments'. Note that here the only real treatment was the experimental creation of pairs in the field. 
# The 'established group' refers to established groups that were captured and recaptured over a similar time interval to the newly created pairs, creating a time-matched comparison. 

mean(experimental$TimeToNextCap[experimental$Treatment == "EstablishedGroup"])
sd(experimental$TimeToNextCap[experimental$Treatment == "EstablishedGroup"])

mean(experimental$TimeToNextCap[experimental$Treatment == "ExperimentalPairing"])
sd(experimental$TimeToNextCap[experimental$Treatment == "ExperimentalPairing"])


t.test(experimental$TimeToNextCap[experimental$Treatment == "EstablishedGroup"], 
       experimental$TimeToNextCap[experimental$Treatment == "ExperimentalPairing"])
# Welch's t=test, t = -1.4674, df = 14.14, p-value = 0.1642

t.test(experimental$WithinGroupRecruits[experimental$Treatment == "EstablishedGroup"], 
       experimental$WithinGroupRecruits[experimental$Treatment == "ExperimentalPairing"])
# Welch's t=test, t = 0.25233, df = 14.022, p-value = 0.8044


#####-----------------  END   ------------------####
