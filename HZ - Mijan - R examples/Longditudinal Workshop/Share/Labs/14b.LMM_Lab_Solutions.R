## Linear Mixed models
## 
## Explore parametric models
## Differences on the random effect assumptions
##
## Use "HIVmetabolites.csv" dataset
## 

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


#Read the data in
hiv.met <- read.csv("https://www.dropbox.com/s/3l5svw86z1s2dt6/HIVmetabolites.csv?dl=1")


#Table summarizing the data
head(hiv.met) # display a few observations
tbl_summary(hiv.met)

## Remove the missing values, not necessary
hiv.met.anl <- hiv.met[!is.na(hiv.met$rel_mid2), ]


## Mixed models in R
library(nlme)
library(help=nlme)

## Pseudo-code for the mixed models: function lme() in the package library(nlme)
#   fixed = outcome~time*factor		## FIXED effects
#   random = ~1+time|id	          ## RANDOM effects
##

## Repeated measures code
rel_mid2.csh.cat <- gls(rel_mid2~factor(time.cat)*nadircd4.cat, 
                       correlation=corCompSymm(form= ~1 | patid),
                       weights=varIdent(form= ~1 | factor(time.cat)),
                       data=hiv.met.anl,
                       control=glsControl(maxIter = 200, msMaxIter = 200, opt="optim"))
summary(rel_mid2.csh.cat)


## Random intercept
fit.RandInt <- lme(rel_mid2 ~ t*nadircd4.cat, # "*" main effects + interaction
                   random=~1|patid,
                   data=hiv.met.anl)
summary(fit.RandInt)


## Random intercept and slope
fit.RandSlp <- lme(rel_mid2 ~ t*nadircd4.cat, # "*" main effects + interaction
                   random=~1+t|patid,
                   data=hiv.met.anl,
                   control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandSlp)

anova(fit.RandSlp, fit.RandInt)

VarCorr(fit.RandSlp)
getVarCov(fit.RandSlp)

## Random intercept and slope
fit.RandSlpQuad <- lme(rel_mid2 ~ t*nadircd4.cat + I(t^2)*nadircd4.cat, 
                   random=~1+t|patid,
                   data=hiv.met.anl,
                   control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandSlpQuad)

# Exercise:
# 1. Fit LMMs for the outcome "abs_cen1" with the "time" and "nadircd4.cat" as predictors
#    a. Random intercepts
#    b. Random intercepts and slopes
# 2. Compare the models from part 1
# 3. Add other baseline predictors to the model
#   a. affecting the starting point (intercepts)
#   b. affecting the trends over time (slopes)

## 1. Random intercept and slope

fit.RandInt.absCen1 <- lme(abs_cen1 ~ t*nadircd4.cat, 
                           random=~1|patid,
                           data=hiv.met.anl, na.action = na.omit,
                           control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandInt.absCen1)

fit.RandSlp.absCen1 <- lme(abs_cen1 ~ t*nadircd4.cat, 
                       random=~1+t|patid,
                       data=hiv.met.anl, na.action = na.omit,
                       control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandSlp.absCen1)

## 2. Compare models
anova(fit.RandInt.absCen1, fit.RandSlp.absCen1)

## 3. Additional covariates
fit.RandSlpBsl.absCen1 <- lme(abs_cen1 ~ t*nadircd4.cat + sexbin + racebin, 
                           random=~1+t|patid,
                           data=hiv.met.anl, na.action = na.omit,
                           control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandSlpBsl.absCen1)


fit.RandSlpTime.absCen1 <- lme(abs_cen1 ~ t*(nadircd4.cat + sexbin + racebin), 
                              random=~1+t|patid,
                              data=hiv.met.anl, na.action = na.omit,
                              control = lmeControl(maxIter = 200, opt="optim"))
summary(fit.RandSlpTime.absCen1)


## Newer library which will be also used in the GLMM portion of the course
library(lme4)

fit4.RandSlp <- lmer(rel_mid2 ~ t*nadircd4.cat + (1+t|patid),
                   data=hiv.met.anl)
summary(fit4.RandSlp) # Note: no p-values
anova(fit4.RandSlp)


## Testing for significance of the coefficients obtained from "lmer"
library(lmerTest)
fit4.RandSlp <- lmer(rel_mid2 ~ t*nadircd4.cat + (1+t|patid),
                     data=hiv.met.anl)
summary(fit4.RandSlp)
anova(fit4.RandSlp)



