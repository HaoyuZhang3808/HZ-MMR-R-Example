##############################
##  2.4 GAMM lab
##############################

################################################################################
# Dataset: "HIVmetabolites.csv"
# Description: Neurocognitive decline of chronically infected HIV patients
#   patid -    patient id
#   sexbin -   male == 1
#   racebin -  white == 1
#   educbin -  educ>12 == 1
#   age     -  age in years
#   cd4val  -  cd4 count at enrollment
#   nadircd4 - lifetime nadir CD4 count
#   durationHIV - duration of HIV infection (years)
#   HAART_years - years on HAART
#   cpescore    - CNS meds penetration score
#   placopy.val - copy number of HIV RNA
#   ADCstage.fac2b - impairment status (impaired == 1)
#   GDS       - Global Deficit Score
#   t           - years since enrollment (continuous)  
#   time.cat - years since enrollment (categorical)
#   abs_*   - absolute metabolite levels in 3 brain regions
#   rel_*   - relative (to Cr) metabolite levels in 3 brain regions
#   rel_mid2 - Mid-frontal Cho/Cr ratio
################################################################################  

library(tidyr)
library(jtools)
library(gtsummary)  # for tables
library(ggplot2)    # for plotting
library(geepack)
library(emmeans)
library(nlme)       # mixed models package
library(HRW)        # datasets 
library(mgcv)       # semiparametric mixed models


data(femSBMD) # Bone mineral density dataset

## Plot by race/ethnicity
ggplot(data = femSBMD,  
       aes(x=age, y=spnbmd, group=as.factor(idnum), color=as.factor(idnum))) +
  geom_line(show.legend=F, alpha=1) +                        #transparency to lines
  scale_x_continuous(breaks=seq(10,25,5)) +          #x-axis 1,2...,6
  labs(title="Bone Mineral Density",
       x = "age (years)",
       y = "spinal bone mineral density (g/cm2)") +  
  facet_wrap(vars(factor(ethnicity))) +      # by category
  geom_smooth(method = loess,                  #adds a smooth line
              aes(group="none"),               #otherwise plots 1 line per ID
              se=F,                            #no 95% confidence band
              color="black") +
  theme_bw()                                   #theme


## Exercise: Plot the outcome "rel_mid2" from the "HIVmetabolites.csv" dataset as as function 
## of age by sex mimicking the code above.

#Read the data
hiv.met <- read.csv("https://www.dropbox.com/s/n96xpn1a25alfcd/HIVmetabolites.csv?dl=1")

summary(hiv.met$age)
hiv.met$age.st <- hiv.met$age + hiv.met$t

## Plot 
ggplot(data = hiv.met,  
       aes(x=age.st, y=rel_mid2, group=as.factor(patid), color=as.factor(patid))) +
  geom_line(show.legend=F, alpha=1) +                        #transparency to lines
  scale_x_continuous(breaks=seq(20,75,5)) +          #x-axis: 20, 25, 30, ..., 70
  labs(title="",
       x = "age (years)",
       y = "AU") +  
  facet_wrap(vars(factor(sexbin, label=c("F","M")))) +      # by category
  geom_smooth(method = loess,                  #adds a smooth line
              aes(group="none"),               #otherwise plots 1 line per ID
              se=F,                            #no 95% confidence band
              color="black") +
  theme_bw()                                   #theme


## Fit a semiparametric model to the "bone mineral density" outcome

fit <- gamm(spnbmd ~ s(age) + black + hispanic + white,
            random = list(idnum = ~1+age), 
            data = femSBMD)

plot(fit$gam, shade = TRUE, shade.col = "palegreen", seWithMean = T)
summary(fit$gam)
intervals(fit$lme)    # confidence intervals for the fixed effects

## Exercise: Build a model for "rel_mid2" with smooth effect of "age" and 
##           parametric effect of "sexbin", "racebin" and "ADCstage.fac2b". 
##           Are there any significant differences in the "Cho/Cr" level?

gamm.hiv <- gamm(rel_mid2 ~ s(age) + racebin + sexbin + ADCstage.fac2b,
            random = list(patid = ~1), 
            data = hiv.met)

plot(gamm.hiv$gam, shade = TRUE, shade.col = "palegreen", seWithMean = T)
summary(gamm.hiv$gam)
intervals(gamm.hiv$lme, which="fixed")    # confidence intervals for the fixed effects





## For illustration purposes only

data(growthIndiana)
growthINblackMales <-
  growthIndiana[(growthIndiana$male == 1)
                & (growthIndiana$black == 1),]

## Plot 
ggplot(data = growthIndiana,  
       aes(x=age, y=height, group=as.factor(idnum), color=as.factor(idnum))) +
  geom_line(show.legend=F, alpha=1) +                        #transparency to lines
  scale_x_continuous(breaks=seq(5,22,3)) +          #x-axis: 20, 25, 30, ..., 70
  labs(title="",
       x = "age (years)",
       y = "height (cm)") +  
  facet_wrap(vars(factor(male, label=c("F","M")))) +      # by category
  geom_smooth(method = loess,                  #adds a smooth line
              aes(group="none"),               #otherwise plots 1 line per ID
              se=F,                            #no 95% confidence band
              color="black") +
  theme_bw()                                   #theme


age <- growthINblackMales$age
height <- growthINblackMales$height
idnum <- growthINblackMales$idnum

uqID <- unique(idnum)
uqID.tab <- table(idnum)
uqID.len <- length(uqID)
growthINblackMales$idnumBM <- as.numeric(
  factor(rep(uqID, uqID.tab),
         labels= 1:uqID.len))
idnumBM <- growthINblackMales$idnumBM


attach(growthINblackMales)
  table(idnum, idnumBM)
detach(growthINblackMales)

  numObs <- length(height)
  numGrp <- uqID.len
  numIntKnotsGbl <- 20
  intKnotsGbl <- quantile(unique(age),
                          seq(0,1,length=numIntKnotsGbl+2))[-c(1,numIntKnotsGbl+2)]
  range.age <- c(5.5,20)
  Zgbl <- ZOSull(age,range.x=range.age,intKnots=intKnotsGbl)
  numIntKnotsGrp <- 10
  intKnotsGrp <- quantile(unique(age),
                          seq(0,1,length=numIntKnotsGrp+2))[-c(1,numIntKnotsGrp+2)]
  Zgrp <- ZOSull(age,range.x=range.age,intKnots=intKnotsGrp)
  
  dummyId <- factor(rep(1,numObs))
  Zblock <- list(dummyId=pdIdent(~-1+Zgbl),
                 idnumBM=pdSymm(~age),
                 idnumBM=pdIdent(~-1+Zgrp)) 
  
  
 
  blkMalGD <- groupedData(height ~ age|rep(1,length = numObs),
                          data = data.frame(height,age,Zgbl,Zgrp,idnumBM))
  fit <- lme(height ~ age,data = blkMalGD,random = Zblock)
  
  
  ng <- 201
  ageg <- seq(range.age[1],range.age[2],length = ng)
  Xg <- cbind(rep(1,ng),ageg)
  Zgblg <- ZOSull(ageg,range.x = range.age,
                  intKnots = intKnotsGbl)
  Zgrpg <- ZOSull(ageg,range.x = range.age,
                  intKnots = intKnotsGrp)
  
  betaHat <- as.vector(fit$coef$fixed)
  uHat <- as.vector(fit$coef$random[[1]])
  fHatg <- as.vector(Xg%*%betaHat + Zgblg%*%uHat)
#  iii. Estimate the subject-specific curves
  curvEsts <- vector("list",numGrp)
  for (i in 1:numGrp)
  {
    uLinHati <- as.vector(fit$coef$random[[2]][i,])
    uSplHati <- as.vector(fit$coef$random[[3]][i,])
    ghati <- Xg%*%uLinHati + Zgrpg%*%uSplHati
    curvEsts[[i]] <- fHatg + ghati
  }

  
## Plot using package "lattice"
library(lattice)
  xyplot(height ~ age|idnumBM,groups = idnumBM,
         data = growthINblackMales,
         strip = FALSE,
         xlab = "age (years)",
         ylab = "height (centimeters)",
         as.table = TRUE, layout = c(4,7),
         panel = function(x,y,subscripts,groups)
         { panel.grid()
           adolNum <- idnumBM[subscripts][1]
           panel.superpose(x,y,subscripts,groups,
                           col = "dodgerblue",type = "b")
           panel.xyplot(ageg,curvEsts[[adolNum]],
                        col = "blue",type = "l")
         })
  
plot(fit)    
summary(fit)

