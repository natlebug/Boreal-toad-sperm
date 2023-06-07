# Hormonal Induction of sperm in Boreal toads
# Script adapted to N. Calatayud's data: Nov-22 by R. Upton
# V1 Nov-22 R. Upton

library(lme4)
library(emmeans)
library(glmmTMB)
library(ggplot2)
library(tidyverse)
library(gridExtra) # Multiple ggplots per page
library(ggeffects)
install.packages('TMB', type = 'source')
library(lmerTest)
library(readxl)
library(bbmle) ## for AICtab
library(dplyr)
library(stats)

# From Ben Bolker's glmm faq page - for assessing overdispersion using Pearson residuals
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

setwd("C:\\Users\\Rose Upton\\OneDrive - The University Of Newcastle\\Documents\\UoN\\Statistics\\Natalie_Calatayud\\Boreal Toad\\Nov-22")
data1=read_excel("RBoreal-2022_dec.xlsx",sheet="hormone treatments")
head(data1)

data1$month_f=factor(data1$month)
data1$year_f=factor(data1$year)
data1$time_f=factor(data1$time)
data1$age_f=factor(data1$age)

# Explore data structure with tables
table(data1$nasrfid)
xtabs(~hormone + nasrfid,data=data1)
table(data1$hormone)
xtabs(~hormone + age,data=data1)
table(data1$season)
table(data1$date)
table(data1$pop)
table(data1$county)

#need to fix different names within spreadsheet if you use this one.
table(data1$source)

#names need fixing
table(data1$caught)

#names need fixing change July to Jul and Jume to Jun so all are consistent - DONE
table(data1$month)

#============================================================================
# Calculate number of sperm/ml
data1 = data1 %>% mutate(sperm_ml = (sperm*dil*10000/pr_sq))

# Scaling factor that converts sperm count to sperm per ml
data1 = mutate(data1,scale_ml=(dil*10000/pr_sq))

# Also reciprocal of the scaling factor per ml and log of reciprocal
data1 = data1 %>% mutate(scale_ml_inv = 1/scale_ml) %>%
  mutate(log_scale_ml_inv = log(scale_ml_inv))
names(data1)

#============================================================================
# Filter information - Remove timepoint 0

# Remove control - all zeroes
data2=filter(data1, time !="0")
nrow(data1)
nrow(data2)

#============================================================================
# Explore the data
# Raw sperm counts
ggplot(data=data2, aes(x=time, y=sperm_ml, color=hormone))+geom_point()
ggplot(data=data2, aes(x=time, y=sperm_ml, color=hormone))+geom_point()+scale_x_log10() 

ggplot(data=data2, aes(x=time, y=sperm_ml, color=hormone))+geom_point()+scale_x_log10() +scale_y_log10()

#Check R calculated sperm/mL = excel calculated sperm/mL - it does, calculation is correct.
ggplot(data=data2, aes(x=sp_mL,y=sperm_ml))+geom_point()

#============================================================================
# Try zero-inflated negative binomial modeling sequence- sperm per ml
# (glm) P -> ZIP -> NB -> ZINB
# P=poisson, ZIP= zero inflated poisson, NB= Negative Binomial, ZINB= Zero inflated negative binomial

#Poisson models
# Runs fine - overdispersed
P_ml <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|nasrfid), data=data2, ziformula=~0, family=poisson)
summary(P_ml)
overdisp_fun((P_ml))

# Runs fine
ZIP1 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|nasrfid), data=data2, ziformula=~1, family=poisson)
summary(ZIP1)

# hour_f to explain zeroes - runs fine
ZIP2 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2)

# hormone to explain zeroes - runs fine
ZIP3 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3)

# year_f to explain zeroes - runs fine
ZIP4 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4)

# month_f to explain zeroes - runs fine
ZIP5 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~month, family=poisson)
summary(ZIP5)

# season to explain zeroes - runs fine
ZIP6 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~season, family=poisson)
summary(ZIP6)

# age to explain zeroes - runs fine
ZIP7 <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7)

# with season as random variable
# Runs fine
ZIP1b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|season), data=data2, ziformula=~1, family=poisson)
summary(ZIP1b)

# hour_f to explain zeroes - runs fine
ZIP2b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2b)

# hormone to explain zeroes - runs fine
ZIP3b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3b)

# year_f to explain zeroes - runs fine
ZIP4b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4b)

# month_f to explain zeroes - runs fine
ZIP5b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~month, family=poisson)
summary(ZIP5b)

# season to explain zeroes - runs fine
ZIP6b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~season, family=poisson)
summary(ZIP6b)

# age to explain zeroes - runs fine
ZIP7b <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7b)

# with month as random variable
# Runs fine
ZIP1c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|month_f), data=data2, ziformula=~1, family=poisson)
summary(ZIP1c)

# hour_f to explain zeroes - runs fine
ZIP2c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2c)

# hormone to explain zeroes - runs fine
ZIP3c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3c)

# year_f to explain zeroes - runs fine
ZIP4c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4c)

# month_f to explain zeroes - runs fine
ZIP5c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~month, family=poisson)
summary(ZIP5c)

# season to explain zeroes - runs fine
ZIP6c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~season, family=poisson)
summary(ZIP6c)

# age to explain zeroes - runs fine
ZIP7c <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7c)


# with year as random variable
# Runs fine
ZIP1d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|year_f), data=data2, ziformula=~1, family=poisson)
summary(ZIP1d)

# hour_f to explain zeroes - runs fine
ZIP2d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2d)

# hormone to explain zeroes - runs fine
ZIP3d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3d)

# year_f to explain zeroes - runs fine
ZIP4d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4d)

# month_f to explain zeroes - runs fine
ZIP5d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~month, family=poisson)
summary(ZIP5d)

# season to explain zeroes - runs fine
ZIP6d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~season, family=poisson)
summary(ZIP6d)

# age to explain zeroes - runs fine
ZIP7d <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7d)

# with age as random variable
# Runs fine
ZIP1e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|age_f), data=data2, ziformula=~1, family=poisson)
summary(ZIP1e)

# hour_f to explain zeroes - runs fine
ZIP2e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2e)

# hormone to explain zeroes - runs fine
ZIP3e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3e)

# year_f to explain zeroes - runs fine
ZIP4e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4e)

# month_f to explain zeroes - runs fine
ZIP5e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~month, family=poisson)
summary(ZIP5e)

# season to explain zeroes - runs fine
ZIP6e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~season, family=poisson)
summary(ZIP6e)

# age to explain zeroes - runs fine
ZIP7e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7e)

# with population as random variable
# Runs fine
ZIP1f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv)  + (1|pop), data=data2, ziformula=~1, family=poisson)
summary(ZIP1f)

# hour_f to explain zeroes - runs fine
ZIP2f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~time_f, family=poisson)
summary(ZIP2f)

# hormone to explain zeroes - runs fine
ZIP3f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~hormone, family=poisson)
summary(ZIP3f)

# year_f to explain zeroes - runs fine
ZIP4f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~year_f, family=poisson)
summary(ZIP4f)

# month_f to explain zeroes - runs fine
ZIP5f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~month, family=poisson)
summary(ZIP5f)

# season to explain zeroes - runs fine
ZIP6f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~season, family=poisson)
summary(ZIP6f)

# age to explain zeroes - runs fine
ZIP7f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~age_f, family=poisson)
summary(ZIP7f)

#check fit of poisson models in comparison to each other:
AICtab(P_ml,ZIP1,ZIP1b,ZIP1c,ZIP1d,ZIP1e,ZIP1f,ZIP2,ZIP2b,ZIP2c,ZIP2d,ZIP2e,ZIP2f,ZIP3,ZIP3b,ZIP3c,ZIP3d,ZIP3e,ZIP3f,
       ZIP4,ZIP4b,ZIP4c,ZIP4d,ZIP4e,ZIP4f,ZIP5,ZIP5b,ZIP5c,ZIP5d,ZIP5e,ZIP5f,ZIP6,ZIP6b,ZIP6c,ZIP6d,ZIP6e,ZIP6f,
       ZIP7,ZIP7b,ZIP7c,ZIP7d,ZIP7e,ZIP7f)


#ZIP models with frog ID have lowest AIC (however residuals are clumped)
simulationOutput <- simulateResiduals(fittedModel = ZIP1, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7, plot = T ,n=1000)

#ZIP models with age as random effect have next best AIC
#residuals slightly clumped and lots of outliers detected
simulationOutput <- simulateResiduals(fittedModel = ZIP1e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7e, plot = T ,n=1000)

#ZIP models with population as random effect have next best fit
#residuals clumped
simulationOutput <- simulateResiduals(fittedModel = ZIP1f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7f, plot = T ,n=1000)

#ZIP models with month as random effect have next best fit
#residuals clumped
simulationOutput <- simulateResiduals(fittedModel = ZIP1c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7c, plot = T ,n=1000)


#season - These residuals are terrible and fail all tests
simulationOutput <- simulateResiduals(fittedModel = ZIP1b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7b, plot = T ,n=1000)

#p_ml model residuals - not too bad but still clumped
simulationOutput <- simulateResiduals(fittedModel = P_ml, plot = T ,n=1000)

#year - These residuals are terrible and fail all tests
#terrible residuals
simulationOutput <- simulateResiduals(fittedModel = ZIP1d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP2d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP3d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP4d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP5d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP6d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZIP7d, plot = T ,n=1000)



# Zero inflated binomials are not sufficient: try series of Negative binomials

#using nasrfid as random variable
# Negative binomial - runs fine
NB  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~0, family=nbinom2)
summary(NB)

# Zero-inflated negative binomial - runs fine
ZINB1  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1)

# hour_f to explain zeroes - runs fine
ZINB2  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2)
overdisp_fun(ZINB2)

# trt to explain zeroes - runs fine
ZINB3  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3)

# year_f to explain zeroes - runs fine
ZINB4  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4)

# month_f to explain zeroes - runs fine
ZINB5  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5)

# season to explain zeroes - runs fine
ZINB6  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6)

# age to explain zeroes - runs fine
ZINB7  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|nasrfid), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7)


#using season as random variable
# Zero-inflated negative binomial - runs fine
ZINB1b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1b)

# hour_f to explain zeroes - runs fine
ZINB2b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2b)

# trt to explain zeroes - runs fine
ZINB3b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3b)

# year_f to explain zeroes - runs fine
ZINB4b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4b)

# month_f to explain zeroes - runs fine
ZINB5b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5b)

# season to explain zeroes - runs fine
ZINB6b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6b)

# age to explain zeroes - runs fine
ZINB7b  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7b)


#using month as random variable
# Zero-inflated negative binomial - runs fine
ZINB1c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1c)

# hour_f to explain zeroes - runs fine
ZINB2c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2c)

# trt to explain zeroes - runs fine
ZINB3c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3c)

# year_f to explain zeroes - runs fine
ZINB4c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4c)

# month_f to explain zeroes - runs fine
ZINB5c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5c)

# season to explain zeroes - runs fine
ZINB6c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6c)

# age to explain zeroes - runs fine
ZINB7c  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|month_f), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7c)



#using year as random variable
# Zero-inflated negative binomial - runs fine
ZINB1d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1d)

# hour_f to explain zeroes - runs fine
ZINB2d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2d)

# trt to explain zeroes - runs fine
ZINB3d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3d)

# year_f to explain zeroes - runs fine
ZINB4d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4d)

# month_f to explain zeroes - runs fine
ZINB5d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5d)

# season to explain zeroes - runs fine
ZINB6d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6d)

# age to explain zeroes - runs fine
ZINB7d  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|year_f), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7d)



#using age as random variable
# Zero-inflated negative binomial - runs fine
ZINB1e  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1e)

# hour_f to explain zeroes - runs fine
ZINB2e  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2e)

# trt to explain zeroes - runs fine
ZINB3e  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3e)

# year_f to explain zeroes - runs fine
ZINB4e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4e)

# month_f to explain zeroes - runs fine
ZINB5e  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5e)

# season to explain zeroes - runs fine
ZINB6e <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6e)

# age to explain zeroes - runs fine
ZINB7e  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|age_f), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7e)



#using population as random variable
# Zero-inflated negative binomial - runs fine
ZINB1f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~1, family=nbinom2)
summary(ZINB1f)

# hour_f to explain zeroes - runs fine
ZINB2f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~time_f, family=nbinom2)
summary(ZINB2f)

# trt to explain zeroes - runs fine
ZINB3f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~hormone, family=nbinom2)
summary(ZINB3f)

# year_f to explain zeroes - runs fine
ZINB4f <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~year_f, family=nbinom2)
summary(ZINB4f)

# month_f to explain zeroes - runs fine
ZINB5f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~month, family=nbinom2)
summary(ZINB5f)

# season to explain zeroes - runs fine
ZINB6f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~season, family=nbinom2)
summary(ZINB6f)

# age to explain zeroes - runs fine
ZINB7f  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|pop), data=data2, ziformula=~age_f, family=nbinom2)
summary(ZINB7f)


#compare the models - different pattern. Here model 3x seems to do best (following those with id as re)
#the 3x models all use hormone to explain zeroes. I think the ZINB models are best for our purposes.
AICtab(NB,ZINB1,ZINB2,ZINB3,ZINB4,ZINB5,ZINB6, ZINB7,ZINB1b,ZINB2b,ZINB3b,ZINB4b,ZINB5b,ZINB6b,ZINB7b,ZINB1c,ZINB2c,ZINB3c,ZINB4c,ZINB5c,ZINB6c,ZINB7c,
       ZINB1d,ZINB2d,ZINB3d,ZINB4d,ZINB5d,ZINB6d,ZINB7d,ZINB1e,ZINB2e,ZINB3e,ZINB4e,ZINB5e,ZINB6e,ZINB7e,ZINB1f,ZINB2f,ZINB3f,ZINB4f,ZINB5f,ZINB6f,ZINB7f)

#ZINB with ID - resdiuals clumped
simulationOutput <- simulateResiduals(fittedModel = ZINB1, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7, plot = T ,n=1000)

#ZINB with age - clumped 7e isn't too bad though
simulationOutput <- simulateResiduals(fittedModel = ZINB1e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6e, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7e, plot = T ,n=1000)

#ZINB with pop - clumped
simulationOutput <- simulateResiduals(fittedModel = ZINB1f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6f, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7f, plot = T ,n=1000)

#ZINB with month - 3c isn't too bad, though has outliers
simulationOutput <- simulateResiduals(fittedModel = ZINB1c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6c, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7c, plot = T ,n=1000)


#season - 3b is quite good, slightly underdispersed, which means it may slightly overestimate
simulationOutput <- simulateResiduals(fittedModel = ZINB1b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6b, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7b, plot = T ,n=1000)

#p_ml model residuals - clumped
simulationOutput <- simulateResiduals(fittedModel = NB, plot = T ,n=1000)

#year - 3d isn't too bad, though it has failed the KS test
simulationOutput <- simulateResiduals(fittedModel = ZINB1d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB2d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB3d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB4d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB5d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB6d, plot = T ,n=1000)
simulationOutput <- simulateResiduals(fittedModel = ZINB7d, plot = T ,n=1000)


#Further investigation of residuals and fit of Model ZINB3b which uses season as random effect and hormone to explain zeroes

#Extract residuals
residuals <- residuals(ZINB3b)

#Plot residuals as a histogram- slightly skewed, but passed the KS test for normality
hist(residuals, main="Histogram of Residuals", xlab="Residuals", col="blue", border="black")

# Extract the residuals and fitted values
residuals <- residuals(ZINB3b)
fitted_values <- fitted(ZINB3b)


# check the length of residuals and fitted_values
# make sure they have the same length
length(residuals) == length(fitted_values)

#line of best fit:
plot(fitted_values,residuals, 
     xlab = "Fitted Values", ylab = "Residuals")
abline(lm(residuals ~ fitted_values))



# Not sure why the predicted function fails (only 171 obs predicted) where the residuals works (178 obs with residuals) 
# but needed a work around to allow plotting the residuals
data2a = select(data2, sperm, sperm_ml,time_f, hormone, log_scale_ml_inv, season)
nrow(data2a)
nrow(na.omit(data2a))

# remove missings so that the predictions can work correctly and residual plots can be added.
data2a = na.omit(data2a)
nrow(data2a)


# Re fit using data3
ZINB3ba  <- glmmTMB(sperm ~ time_f + hormone + offset(log_scale_ml_inv) + (1|season), data=data2a, ziformula=~hormone, family=nbinom2)



# Traditional residual diagnostic plot residuals on y axis against predicted values on the x axis
data2a = mutate(data2a, 
                resid.response = residuals(ZINB3b,type="response"),
                pred.response = predict(ZINB3b,type="response"),
                resid.pearson = residuals(ZINB3b,type="pearson"),
                pred.link = predict(ZINB3b,type="link"),
                resid.working = residuals(ZINB3b,type="working"))

# Check residuals - most important plot I think is the pearson residual (Kim C)
# Pearson residuals (middle plot show data points more than 1 sd above mean)
op=par(mfrow=c(1,3))
plot(resid.response ~ pred.response,data=data2a)
abline(h=0)
plot(resid.pearson ~ pred.link,data=data2a)
abline(h=0)
plot(resid.working ~ pred.link,data=data2a)
abline(h=0)
par(op)

#Check fit of model against raw data- Fritz says this is the biggest test of a model
# This model fits the data pretty well. Due to a slightly underdispersed model, we are slightly overpredicting the means when compared to raw data.


ggplot(data=data2, aes(x=time, y=sperm_ml, color=hormone))+geom_point()+scale_x_log10() +scale_y_log10()
emmip(ZINB3b, hormone ~ time_f  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")

#Underdispersed?
simulationOutput <- simulateResiduals(fittedModel = ZINB3b, plot = T ,n=1000)
testDispersion(simulationOutput,alternative = "less")
#-------------Means from raw data----------------------------------------------
# # Calculate the means of column sp_mL grouped by columns hormone and time_f
# agg_means <- aggregate(sp_mL ~ hormone + time_f, data = data2, mean)

# Group the data by the columns hormone and time_f
grouped_data <- with(data2, split(sp_mL, list(hormone, time_f)))

# Calculate the means and standard deviations of column sp_mL for each group
agg_result <- lapply(grouped_data, function(x) {
  c(mean = mean(x), sd = sd(x))
})
agg_result <- as.data.frame(do.call(rbind, agg_result))
agg_result

#Compare raw means to model estimated marginal means means
rawvmodel=read_excel("anal1_raw_v_model.xlsx",sheet="Sheet1")
ggplot(data=rawvmodel, aes(x=raw.mean,y=model.mean))+geom_point()




#===================================================================================================================================================
# Best model is ZINB3b - time and hormone is highly significant
summary(ZINB3b)
#this function is not ideal for ZINB models, but shows good dispersion.
overdisp_fun(ZINB3b)
drop1(ZINB3b,test="Chisq")

# # Plot count and zeros models - version 1
count2=emmip(ZINB3b, ~ time_f | hormone, CIs=T,type="response",component="cond",offset=0) + labs(title="Counts adjusted for zero inflation",x="Time",y="Counts")
# # Proportion of Zero inflation
zero2=emmip(ZINB3b.g, ~hormone, CIs=T,type="response",component="zi",offset=0) + labs(title="Proportion of excess zeros",x="Month",y="Proportion")
grid.arrange(count2, zero2, ncol=2)

#Two graphs for paper
emmip(ZINB3b, hormone ~ time_f  ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")

emmip(ZINB3b, ~ time_f | hormone ,CIs=T,type="response",component="cond",offset=0) + labs(x="Hours",y="Sperm Concentration (cells/mL)")

# Estimated Marginal Means and 95% CIs
ZINB3b.g=ref_grid(ZINB3b)
emmeans(ZINB3b.g,"hormone",by="time_f", type="response", offset=0)
emmeans(ZINB3b.g,"time_f",by="hormone", type="response", offset=0)


#Pairwise comparison of ratios
# hormone
emmeans(ZINB3b.g,pairwise~hormone,type="link")
pairs.hormone=emmeans(ZINB3b.g,pairwise~hormone,type="response",adjust="none")
confint(pairs.hormone)

emmeans(ZINB3b.g,revpairwise~hormone,type="link")
pairs.hormone=emmeans(ZINB3b.g,revpairwise~hormone,type="response",adjust="none")
confint(pairs.hormone)

# Hour post injection (time_f)
emmeans(ZINB3b.g,pairwise~time_f,type="link")
pairs.time_f=emmeans(ZINB3b.g,pairwise~time_f,type="response",adjust="none")
confint(pairs.time_f)

emmeans(ZINB3b.g,revpairwise~time_f,type="link")
pairs.time_f=emmeans(ZINB3b.g,revpairwise~time_f,type="response",adjust="none")
confint(pairs.time_f)

#-------------------------------------------------------------------------------
