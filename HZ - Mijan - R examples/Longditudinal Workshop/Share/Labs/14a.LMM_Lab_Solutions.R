# --------------------------------------------------------------------------
#  Linear mixed models LAB
# --------------------------------------------------------------------------

    # Dataset: Strength 
    #
    # Study: Subjects in an exercise therapy study were assigned to one of two 
    #        weightlifting programs. 
    #        First program (treatment 1), the number of repetitions was increased 
    #        Second program (treatment 2), the weight was increased
    #        For illustration, we focus only on measures of strength taken at:
    #        baseline (day 0) and on days 4, 6, 8, and 12.
    #
    # Variables:
    #    y - strength measurement
    #    id - subject's id
    #    time - time in days
    #    trt - exercise group



library(tidyverse)
library(gtsummary)   #for tables
library(tidyr)       #for the function pivot_longer and pivot_wider
library(ggplot2)     #for plots
library(nlme)        #for fitting longitudinal models
library(GGally)      #for pairs plot


## Read the data in
stren <- read.table("https://www.dropbox.com/s/f1n20of73hrispo/stren.txt?dl=1", 
                      header=F, na.strings = ".")

names(stren) <- c("id", "trt", paste("y", 1:7, sep=""))

## Print the first few observations
head(stren)

## Pairs plot: all longitudinal observations
ggpairs(stren[,3:9])


stren.mean <- apply(stren[,3:9], 2, mean, na.rm=T)
matplot(stren.mean, type='l')

stren.uni <- data.frame(id=rep(stren$id, each=7),
                      trt=rep(stren$trt, each=7),
                      y=as.numeric(t(as.matrix(stren[,3:9]))),
                      time=rep(seq(0,12,2)),
                      time.cat=rep(1:7))
stren.uni[1:14,]

# Create an unevenly-spaced dataset
stren.uni.sel <- stren.uni[stren.uni$time %in% c(0,4,6,8,12), ]
stren.uni.sel <- stren.uni.sel[!is.na(stren.uni.sel$y), ]

tjj.mean <- tapply(stren.uni.sel$y, list(stren.uni.sel$trt, stren.uni.sel$time), mean)
matplot(c(0,4,6,8,12), t(tjj.mean), type='b', ylim=c(78,85),
        xlab="Time", ylab="Strength")

## Spaghetti plot
ggplot(data = stren.uni.sel,  
       aes(x=time, y=y, group=id, color=factor(id))) +
  scale_x_continuous(breaks=2*c(0:6)) +          #x-axis 0,2,4,...,12
  facet_grid(cols=vars(trt)) +
  stat_summary(fun=mean, geom="line",show.legend = F, lwd=3, aes(group="none")) +
  geom_line(show.legend = F, alpha=.5) 

## Mean plot
ggplot(stren.uni.sel, aes(x=time, y=y, group=factor(trt), color=factor(trt))) + 
  scale_x_continuous(breaks=2*c(0:6)) +          #x-axis 0,2,4,...,12
  stat_summary(fun=mean, geom="line") +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar") + 
  coord_cartesian(ylim=c(76,86))


## Unstructured covariance matrix
stren.un <- gls(y~factor(time)*trt, 
                correlation=corSymm(form= ~1 | factor(id)),
                weights=varIdent(form= ~1 | factor(time)),
                data=stren.uni.sel)
summary(stren.un)

## Estimated covariance and correlation matrices
stren.un.VarCov <- getVarCov(stren.un)
stren.un.VarCov

sqrt.un.VarCov <- sqrt(outer(diag(stren.un.VarCov), diag(stren.un.VarCov)))
stren.un.Corr <- stren.un.VarCov/sqrt.un.VarCov
stren.un.Corr

## AR(1) discrete time model
stren.ar1 <- gls(y~factor(time)*trt, 
                correlation=corAR1(form= ~1 | factor(id)),
                data=stren.uni.sel)
summary(stren.ar1)

stren.ar1.VarCov <- getVarCov(stren.ar1)
stren.ar1.VarCov
sqrt.ar1.VarCov <- sqrt(outer(diag(stren.ar1.VarCov), diag(stren.ar1.VarCov)))
stren.ar1.Corr <- stren.ar1.VarCov/sqrt.ar1.VarCov
stren.ar1.Corr


## AR(1) continuous time model
stren.car1 <- gls(y~factor(time)*trt, 
                 correlation=corCAR1(form= ~time | factor(id)),
                 data=stren.uni.sel)
summary(stren.car1)

stren.car1.VarCov <- getVarCov(stren.car1)
stren.car1.VarCov
sqrt.car1.VarCov <- sqrt(outer(diag(stren.car1.VarCov), diag(stren.car1.VarCov)))
stren.car1.Corr <- stren.car1.VarCov/sqrt.car1.VarCov
stren.car1.Corr


## Mean model - ML estimation
stren.car1.mean.ml <- gls(y~factor(time)*trt, 
                correlation=corCAR1(form= ~time | factor(id)),
                data=stren.uni.sel, method="ML")
summary(stren.car1.mean.ml)

stren.car1.lin.ml <- gls(y~time*trt, 
                  correlation=corCAR1(form= ~time | factor(id)),
                  data=stren.uni.sel, method="ML")
summary(stren.car1.lin.ml)

anova(stren.car1.mean.ml, stren.car1.lin.ml)


## Mean model - REML
stren.car1.lin.reml <- gls(y~time*trt, 
                         correlation=corCAR1(form= ~time | factor(id)),
                         data=stren.uni.sel, method="REML")
summary(stren.car1.lin.reml)


stren.car1.lin.reml.trt <- gls(y~ time + time:trt, 
                           correlation=corCAR1(form= ~time | factor(id)),
                           data=stren.uni.sel, method="REML")
summary(stren.car1.lin.reml.trt)


## Slide 14: Random intercept and random slope model
stren.lme.int.slp <- lme(y~time*trt, 
                         random=~1+time|factor(id),
                         data=stren.uni.sel)
summary(stren.lme.int.slp)

getVarCov(stren.lme.int.slp)
cov2cor(simplify2array(getVarCov(stren.lme.int.slp,type="marginal")[[1]]))


## Slide 19: (alternatives to lme() --> gls() )
stren.lme.int <- lme(y~time*trt, 
                     random=~1|factor(id),
                     data=stren.uni.sel)
summary(stren.lme.int)

stren.cs <- gls(y~time*trt, 
                correlation=corCompSymm(form= ~time | factor(id)),
                data=stren.uni.sel)
summary(stren.cs)


## Slide 25:
coef(stren.lme.int.slp)
stren.lme.int.slp$coef$random
stren.lme.int.slp$coef$fixed

## Slide 28: fitted values
stren.lme.int.slp.pred <- cbind(stren.uni.sel, 
                                predict(stren.lme.int.slp),
                                predict(stren.lme.int.slp, level = 0))
names(stren.lme.int.slp.pred)[6:7] <- c("pred.int.slp", "pred.mean")

head(stren.lme.int.slp.pred)

ggplot(data = stren.lme.int.slp.pred,  
       aes(x=time, y=pred.int.slp, group=id, color=factor(id))) +
  scale_x_continuous(breaks=2*c(0:6)) +          #x-axis 0,2,4,...,12
  geom_line(show.legend = F, alpha=.5) +
  facet_grid(cols=vars(trt)) + 
  geom_line(aes(y=pred.mean), lwd=2, show.legend=F)

##
## Exercise: Reproduce the plot of the fitted lines for the "random intercept"-only model
##  

## Fitted values: random intercept model
stren.lme.int.pred <- cbind(stren.uni.sel, 
                                predict(stren.lme.int),
                                predict(stren.lme.int, level = 0))
names(stren.lme.int.pred)[6:7] <- c("pred.int", "pred.mean")

head(stren.lme.int.pred)

ggplot(data = stren.lme.int.pred,  
       aes(x=time, y=pred.int, group=id, color=factor(id))) +
  scale_x_continuous(breaks=2*c(0:6)) +          #x-axis 0,2,4,...,12
  geom_line(show.legend = F, alpha=.5) +
  facet_grid(cols=vars(trt)) + 
  geom_line(aes(y=pred.mean), lwd=2, show.legend=F)

