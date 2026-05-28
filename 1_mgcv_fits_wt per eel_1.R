##--------------------------------
## Fit trend and covariate models
## NK: Feb 27, 2026, updated from CM
##--------------------------------

########  year is continuous smoothed (not a factor)  ###################

########  In this code, we are standardizing the year effect for temperature  ###################
########  _yrstempd means yr as smooth, watertemp as density covariate  ###################

##--------------------------------
## Fit trend and covariate models
##--------------------------------
library(mgcv)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")


## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

#### WITH environmental covariates, years as smooth ####

##-----------------
## WEIGHT DATA FITS 
##-----------------
#containers
weight_fits_yrstempd <- list()
weight_effects_pred_yrstempd <- NULL
weight_year_pred_yrstempd <- NULL

for(lake in lakes){
  print(lake)
  sub_dat <- subset(weightenvdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake %in% c("BOH", "Bunaveela")){
    form <- as.formula(wt ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other))
                         #guard + 
                         #offset(log(Effort)))
  }
  if(lake == "Furnace"){
    form <- as.formula(wt ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other)) 
                         #guard + 
                         #offset(log(Effort)))
  }
  if(lake == "Feeagh"){
    form <- as.formula(wt ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other)) 
                         #guard + 
                         #offset(log(Effort)))
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  #gam.vcomp(f0)
  print(AIC(f0))
  print(summary(f0))
  ###
  weight_fits_yrstempd[[lake]] <- f0
  ## effect predictions 
  m <- 100
  pred_df0 <- data.frame(
    Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
    DOY = seq(min(sub_dat$DOY), max(sub_dat$DOY), length = m),
    fSite = unique(sub_dat$fSite)[1],
    watertemp = seq(min(sub_dat$watertemp, na.rm=TRUE), max(sub_dat$watertemp, na.rm=TRUE), length = m),
    pressure = seq(min(sub_dat$pressure, na.rm=TRUE), max(sub_dat$pressure, na.rm=TRUE), length = m),
    wind = seq(min(sub_dat$wind, na.rm=TRUE), max(sub_dat$wind, na.rm=TRUE), length = m),
    clouds = seq(min(sub_dat$clouds, na.rm=TRUE), max(sub_dat$clouds, na.rm=TRUE), length = m),
    moonlight = seq(min(sub_dat$moonlight, na.rm=TRUE), max(sub_dat$moonlight, na.rm=TRUE), length = m),
    waterlev = seq(min(sub_dat$waterlev, na.rm=TRUE), max(sub_dat$waterlev, na.rm=TRUE), length = m))
    #Effort = 1) ## setting this to 1 net-night
  ## predictions for non random effects
  pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
  zs <- colnames(pred0$fit)
  ## for Feeagh only
  ##zs <- zs[zs != "s(trap_number):surveyIFI"]
  zs <- zs[zs != "ti(clouds,moonlight)"]
  lake_effects <- NULL
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    if(v == "trap_number:surveyRussell"){
      v <- "trap_number"
    }
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      x = pred_df0[,v],
                      yhat = pred0$fit[, z],
                      ylwr = pred0$fit[, z] - 2 * pred0$se.fit[, z],
                      yupr = pred0$fit[, z] + 2 * pred0$se.fit[, z]
    )
    tmp <- unique(tmp)
    lake_effects <- rbind(lake_effects, tmp)
    rm(tmp)
  }
  ## predictions for site
  pred_df1 <- expand.grid(
    Year = sub_dat$Year[1],
    DOY = mean(sub_dat$DOY),
    #Effort = 1, ## setting this to one net
    fSite = unique(sub_dat$fSite),
    watertemp = mean(sub_dat$watertemp, na.rm=TRUE),
    pressure = mean(sub_dat$pressure, na.rm=TRUE),
    wind = mean(sub_dat$wind, na.rm=TRUE),
    clouds = mean(sub_dat$clouds, na.rm=TRUE),
    moonlight = mean(sub_dat$moonlight, na.rm=TRUE),
    waterlev = mean(sub_dat$waterlev, na.rm=TRUE))
  ## predictions for random effects
  pred1 <- predict(f0, newdata = pred_df1, type = "terms", se.fit = TRUE)
  zs <- c("s(fSite)")
  random_effects <- NULL
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      x = pred_df1[,v],
                      yhat = pred1$fit[, z],
                      ylwr = pred1$fit[, z] - 2 * pred1$se.fit[, z],
                      yupr = pred1$fit[, z] + 2 * pred1$se.fit[, z]
    )
    tmp <- unique(tmp)
    random_effects <- rbind(random_effects, tmp)
    rm(tmp)
  }
  all_effects <- rbind(lake_effects, random_effects)
  weight_effects_pred_yrstempd <- rbind(weight_effects_pred_yrstempd, all_effects)
  ##------------------------
  ## GET YEARLY PREDICTIONS--now standardized by setting temperature to year fitted values
  ##------------------------
  ## here we set non-year and temperature continuous covariates to their mean and exclude random effects
  # 1) Compute means for all continuous covariates (except temp)
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
      #trap_depth = mean(trap_depth, na.rm = TRUE),
      #trap_gradient = mean(trap_gradient, na.rm = TRUE),
      pressure = mean(pressure, na.rm = TRUE),
      wind = mean(wind, na.rm = TRUE),
      clouds = mean(clouds, na.rm = TRUE),
      moonlight = mean(moonlight, na.rm = TRUE),
      waterlev = mean(waterlev, na.rm = TRUE)
    )
  # 2) Compute YEAR-SPECIFIC mean temperature
  year_temp <- sub_dat %>%
    group_by(Year) %>%
    summarise(watertemp = mean(watertemp, na.rm = TRUE), .groups = "drop")
  
  pred_df2 <- data.frame(Year = year_temp$Year,
                         DOY = means$DOY,
                         fSite = unique(sub_dat$fSite)[1],
                         watertemp = year_temp$watertemp,  # <-- CONDITIONED on actual temps
                         pressure = means$pressure,
                         wind = means$wind,
                         clouds = means$clouds,
                         moonlight = means$moonlight,
                         waterlev = means$waterlev)
                         #Effort = 1)
  
  
  pred2 <- predict(f0, newdata = pred_df2, type = "terms", terms = "s(Year)", se.fit = TRUE,
                   ##exclude = c("s(fSite)", "s(DOY)"))
                   exclude = c("s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$yhat <- exp(pred2$fit)
  pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
  pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
  weight_year_pred_yrstempd <- rbind(weight_year_pred_yrstempd, pred_year)
}

save(weight_year_pred_yrstempd, file = "weight_year_pred_yrstempd.RData")
save(weight_effects_pred_yrstempd, file = "weight_effects_pred_yrstempd.RData")
save(weight_fits_yrstempd, file = "weight_fits_yrstempd.RData")

