##--------------------------------
## Fit trend and covariate models
## NK: Feb 27, 2026, updated from CM
##--------------------------------

########  this file contains:
###### year as smooth, no environmental covariates included (_1)--NOT CURRENTLY USED in plots

###### environmental covariates included, not conditioned on temperature (temp is considered a catchability, not density variable) #######
########  _yrstempc means yr as smooth, watertemp as capacity covariate  ###################

##--------------------------------
## Fit trend and covariate models
##--------------------------------
library(mgcv)
library(ggplot2); theme_set(theme_bw())

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

## models WITH environmental covariates
#### not conditioned on temperature (temp is considered a catchability, not density variable) #######
########  _yrstempc means yr as smooth, watertemp as capacity  ###################

##-----------------
## COUNT DATA FITS 
##-----------------
## containers
count_fits_yrstempc <- list()
effects_pred_yrstempc <- NULL
year_pred_yrstempc <- NULL

for(lake in lakes){
  print(lake)
  ## remove October sampling, which was out of the sampling season 
  sub_dat <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake == "Feeagh"){
    form <- as.formula(count ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         s(trap_depth, m = 2, bs = bs_other) +
                         s(trap_gradient, m = 2, bs = bs_other) +
                         s(trap_number, by = survey, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard +
                         offset(log(Effort)))
  }
  if(lake == "BOH"){
    form <- as.formula(count ~ 
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, k = 5, m = 2, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, m = 2, bs = bs_other) +
                         s(trap_gradient, k = 5, m = 2, bs = bs_other) +
                         s(trap_number, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard +
                         offset(log(Effort)))
  }
  if(lake %in% c("Furnace", "Bunaveela")){
    form <- as.formula(count ~ 
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, k = 5, m = 2, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, k= 5, m = 2, bs = bs_other) +
                         s(trap_gradient, k= 5, m = 2, bs = bs_other) +
                         s(trap_number, k= 5, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         s(watertemp, m = 2, bs = bs_other) +
                         s(pressure, m = 2, bs = bs_other) +
                         s(wind, m = 2, bs = bs_other) +
                         s(clouds, m = 2, bs = bs_other) +
                         s(moonlight, m = 2, bs = bs_other) +
                         ti(clouds, moonlight, m = 2, bs = bs_other) +
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard +
                         offset(log(Effort)))
  }
  ## fit the model
  f0 <- gam(form,
            select = TRUE,
            method = "REML",
            family = nb(),
            data = sub_dat)
  #gam.vcomp(f0)
  ##
  count_fits_yrstempc[[lake]] <- f0
  ## effect predictions
  m <- 100
  pred_df0 <- data.frame(
    Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
    DOY = seq(min(sub_dat$DOY), max(sub_dat$DOY), length = m),
    trap_depth = seq(min(sub_dat$trap_depth), max(sub_dat$trap_depth), length = m),
    trap_gradient = seq(min(sub_dat$trap_gradient), max(sub_dat$trap_gradient), length = m),
    trap_number = rep(1:20, 5),
    survey = "Russell",
    fchain = unique(sub_dat$fchain)[1],
    fSite = unique(sub_dat$fSite)[1],
    watertemp = seq(min(sub_dat$watertemp, na.rm=TRUE), max(sub_dat$watertemp, na.rm=TRUE), length = m),
    pressure = seq(min(sub_dat$pressure, na.rm=TRUE), max(sub_dat$pressure, na.rm=TRUE), length = m),
    wind = seq(min(sub_dat$wind, na.rm=TRUE), max(sub_dat$wind, na.rm=TRUE), length = m),
    clouds = seq(min(sub_dat$clouds, na.rm=TRUE), max(sub_dat$clouds, na.rm=TRUE), length = m),
    moonlight = seq(min(sub_dat$moonlight, na.rm=TRUE), max(sub_dat$moonlight, na.rm=TRUE), length = m),
    waterlev = seq(min(sub_dat$waterlev, na.rm=TRUE), max(sub_dat$waterlev, na.rm=TRUE), length = m),
    Effort = 1)
  if(lake == "Feeagh"){
    ## two survey trap number effects in Feeagh
    tmp <- pred_df0
    tmp <- subset(tmp, trap_number <= 10)
    tmp$survey <- "IFI"
    pred_df0 <- rbind(pred_df0, tmp)
  }    
  ## predictions for non random effects
  pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
  zs <- colnames(pred0$fit)
  ## for Feeagh only
  zs <- zs[zs != "s(trap_number):surveyIFI"]
  zs <- zs[zs != "ti(clouds,moonlight)"]
  lake_effects <- NULL
  idx <- which(pred_df0$survey == "Russell")
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    if(v == "trap_number:surveyRussell"){
      v <- "trap_number"
    }
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "Russell",
                      x = pred_df0[idx,v],
                      yhat = pred0$fit[idx, z],
                      ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                      yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
    )
    tmp <- unique(tmp)
    lake_effects <- rbind(lake_effects, tmp)
    rm(tmp)
  }
  ## add in IFI on Feeagh
  if(lake == "Feeagh"){
    z <- "s(trap_number):surveyIFI"
    idx <- which(pred_df0$survey == "IFI")
    v <- "trap_number"
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "IFI",
                      x = pred_df0[idx,v],
                      yhat = pred0$fit[idx, z],
                      ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                      yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
    )
    tmp <- unique(tmp)
    lake_effects <- rbind(lake_effects, tmp)
  }
  ## predictions for site and chain
  pred_df1 <- expand.grid(
    Year = sub_dat$Year[1],
    DOY = mean(sub_dat$DOY),
    trap_depth = mean(sub_dat$trap_depth),
    trap_gradient = mean(sub_dat$trap_depth),
    trap_number = 10,
    survey = "Russell",
    fchain = unique(sub_dat$fchain),
    fSite = unique(sub_dat$fSite),
    watertemp = mean(sub_dat$watertemp, na.rm=TRUE),
    pressure = mean(sub_dat$pressure, na.rm=TRUE),
    wind = mean(sub_dat$wind, na.rm=TRUE),
    clouds = mean(sub_dat$clouds, na.rm=TRUE),
    moonlight = mean(sub_dat$moonlight, na.rm=TRUE),
    waterlev = mean(sub_dat$waterlev, na.rm=TRUE),
    Effort = 1)
  ## predictions for random effects
  pred1 <- predict(f0, newdata = pred_df1, type = "terms", se.fit = TRUE)
  zs <- c("s(fSite)", "s(fchain)")
  random_effects <- NULL
  ## 
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "Russell",
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
  effects_pred_yrstempc <- rbind(effects_pred_yrstempc, all_effects)
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                         DOY = mean(sub_dat$DOY, na.rm = TRUE),
                         trap_depth = mean(sub_dat$trap_depth, na.rm = TRUE),
                         trap_gradient = mean(sub_dat$trap_gradient, na.rm = TRUE),
                         trap_number = 10,
                         survey = "Russell",
                         fchain = unique(sub_dat$fchain)[1],
                         fSite = unique(sub_dat$fSite)[1],
                         watertemp = mean(sub_dat$watertemp, na.rm = TRUE),
                         pressure = mean(sub_dat$pressure, na.rm = TRUE),
                         wind = mean(sub_dat$wind, na.rm = TRUE),
                         clouds = mean(sub_dat$clouds, na.rm = TRUE),
                         moonlight = mean(sub_dat$moonlight, na.rm = TRUE),
                         waterlev = mean(sub_dat$waterlev, na.rm=TRUE),
                         Effort = 2) ## 2 codends is a fyke net
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   exclude = c("s(fchain)", "s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$yhat <- exp(pred2$fit)
  pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
  pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
  year_pred_yrstempc <- rbind(year_pred_yrstempc, pred_year)
}

save(year_pred_yrstempc, file = "year_pred_yrstempc.RData")
save(effects_pred_yrstempc, file = "effects_pred_yrstempc.RData")
save(count_fits_yrstempc, file = "count_fits_yrstempc.RData")



##-----------------
## WEIGHT DATA FITS 
##-----------------
##MODEL WITH ENVIRONMENTAL COVARIATES
#all variables as fixed effects, no lagging

#containers
weight_fits_yrstempc <- list()
weight_effects_pred_yrstempc <- NULL
weight_year_pred_yrstempc <- NULL

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
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard + 
                         offset(log(Effort))
    )
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
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard + 
                         offset(log(Effort))
    )
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
                         s(waterlev, m = 2, bs = bs_other) +
                         #guard + 
                         offset(log(Effort))
    )
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  #gam.vcomp(f0)
  ###
  weight_fits_yrstempc[[lake]] <- f0
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
    waterlev = seq(min(sub_dat$waterlev, na.rm=TRUE), max(sub_dat$waterlev, na.rm=TRUE), length = m),
    Effort = 1) ## setting this to 1 net-night
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
    Effort = 1, ## setting this to one net
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
  weight_effects_pred_yrstempc <- rbind(weight_effects_pred_yrstempc, all_effects)
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                         DOY = 200,
                         fSite = unique(sub_dat$fSite)[1],
                         watertemp = mean(sub_dat$watertemp, na.rm = TRUE),
                         pressure = mean(sub_dat$pressure, na.rm = TRUE),
                         wind = mean(sub_dat$wind, na.rm = TRUE),
                         clouds = mean(sub_dat$clouds, na.rm = TRUE),
                         moonlight = mean(sub_dat$moonlight, na.rm = TRUE),
                         waterlev = mean(sub_dat$waterlev, na.rm=TRUE),
                         Effort = 1) ## setting this to one net
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   ##exclude = c("s(fSite)", "s(DOY)"))
                   exclude = c("s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$yhat <- exp(pred2$fit)
  pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
  pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
  weight_year_pred_yrstempc <- rbind(weight_year_pred_yrstempc, pred_year)
}

save(weight_year_pred_yrstempc, file = "weight_year_pred_yrstempc.RData")
save(weight_effects_pred_yrstempc, file = "weight_effects_pred_yrstempc.RData")
save(weight_fits_yrstempc, file = "weight_fits_yrstempc.RData")










##MODEL WITH NO ENVIRONMENTAL COVARIATES, year as smooth (_1)

##-----------------
## COUNT DATA FITS 
##-----------------
load("all_lcdat.RData")

## containers
count_fits_1 <- list()
effects_pred_1 <- NULL
year_pred_1 <- NULL

for(lake in lakes){
  print(lake)
  ## remove October sampling, which was out of the sampling season 
  sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake == "Feeagh"){
    form <- as.formula(count ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         s(trap_depth, m = 2, bs = bs_other) +
                         s(trap_gradient, m = 2, bs = bs_other) +
                         s(trap_number, by = survey, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         #guard +
                         offset(log(Effort)))
  }
  if(lake == "BOH"){
    form <- as.formula(count ~ 
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, k = 5, m = 2, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, m = 2, bs = bs_other) +
                         s(trap_gradient, k = 5, m = 2, bs = bs_other) +
                         s(trap_number, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         #guard +
                         offset(log(Effort)))
  }
  if(lake %in% c("Furnace", "Bunaveela")){
    form <- as.formula(count ~ 
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, k = 5, m = 2, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, k= 5, m = 2, bs = bs_other) +
                         s(trap_gradient, k= 5, m = 2, bs = bs_other) +
                         s(trap_number, k= 5, m = 2, bs = bs_other) +
                         s(fchain, bs = "re") +
                         #guard +
                         offset(log(Effort)))
  }
  ## fit the model
  f0 <- gam(form,
            select = TRUE,
            method = "REML",
            family = nb(),
            data = sub_dat)
  #plot(f0)
  #gam.vcomp(f0)
  ###
  count_fits_1[[lake]] <- f0
  ## effect predictions
  m <- 100
  pred_df0 <- data.frame(
    Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
    DOY = seq(min(sub_dat$DOY), max(sub_dat$DOY), length = m),
    trap_depth = seq(min(sub_dat$trap_depth), max(sub_dat$trap_depth), length = m),
    trap_gradient = seq(min(sub_dat$trap_gradient), max(sub_dat$trap_gradient), length = m),
    trap_number = rep(1:20, 5),
    survey = "Russell",
    fchain = unique(sub_dat$fchain)[1],
    fSite = unique(sub_dat$fSite)[1],
    Effort = 1)
  if(lake == "Feeagh"){
    ## two survey trap number effects in Feeagh
    tmp <- pred_df0
    tmp <- subset(tmp, trap_number <= 10)
    tmp$survey <- "IFI"
    pred_df0 <- rbind(pred_df0, tmp)
  }    
  ## predictions for non random effects
  pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
  zs <- colnames(pred0$fit)
  ## for Feeagh only
  zs <- zs[zs != "s(trap_number):surveyIFI"]
  lake_effects <- NULL
  idx <- which(pred_df0$survey == "Russell")
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    if(v == "trap_number:surveyRussell"){
      v <- "trap_number"
    }
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "Russell",
                      x = pred_df0[idx,v],
                      yhat = pred0$fit[idx, z],
                      ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                      yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
    )
    tmp <- unique(tmp)
    lake_effects <- rbind(lake_effects, tmp)
    rm(tmp)
  }
  ## add in IFI on Feeagh
  if(lake == "Feeagh"){
    z <- "s(trap_number):surveyIFI"
    idx <- which(pred_df0$survey == "IFI")
    v <- "trap_number"
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "IFI",
                      x = pred_df0[idx,v],
                      yhat = pred0$fit[idx, z],
                      ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                      yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
    )
    tmp <- unique(tmp)
    lake_effects <- rbind(lake_effects, tmp)
  }
  ## predictions for site and chain
  pred_df1 <- expand.grid(
    Year = sub_dat$Year[1],
    DOY = mean(sub_dat$DOY),
    trap_depth = mean(sub_dat$trap_depth),
    trap_gradient = mean(sub_dat$trap_depth),
    trap_number = 10,
    survey = "Russell",
    fchain = unique(sub_dat$fchain),
    fSite = unique(sub_dat$fSite),
    Effort = 1)
  ## predictions for random effects
  pred1 <- predict(f0, newdata = pred_df1, type = "terms", se.fit = TRUE)
  zs <- c("s(fSite)", "s(fchain)")
  random_effects <- NULL
  ## 
  for(z in zs){
    v <- gsub("(s\\(|\\))", "", z)
    tmp <- data.frame(Lake = lake,
                      variable = v,
                      survey = "Russell",
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
  effects_pred_1 <- rbind(effects_pred_1, all_effects)
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                         DOY = mean(sub_dat$DOY, na.rm = TRUE),
                         trap_depth = mean(sub_dat$trap_depth, na.rm = TRUE),
                         trap_gradient = mean(sub_dat$trap_gradient, na.rm = TRUE),
                         trap_number = 10,
                         survey = "Russell",
                         fchain = unique(sub_dat$fchain)[1],
                         fSite = unique(sub_dat$fSite)[1],
                         Effort = 2) ## 2 codends is a fyke net
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   exclude = c("s(fchain)", "s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$yhat <- exp(pred2$fit)
  pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
  pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
  year_pred_1 <- rbind(year_pred_1, pred_year)
}

save(year_pred_1, file = "year_pred.RData")
save(effects_pred_1, file = "effects_pred.RData")
save(count_fits_1, file = "count_fits.RData")


##-----------------
## WEIGHT DATA FITS 
##-----------------
load("wdat.RData")

#containers
weight_fits_1 <- list()
weight_effects_pred_1 <- NULL
weight_year_pred_1 <- NULL

for(lake in lakes){
  print(lake)
  sub_dat <- subset(wdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake != "Furnace"){
    form <- as.formula(wt ~
                         s(Year, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         #guard + 
                         offset(log(Effort))
    )
  }else{
    form <- as.formula(wt ~
                         s(Year, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") +
                         #guard + 
                         offset(log(Effort))
    )
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  #gam.vcomp(f0)
  ###
  weight_fits_1[[lake]] <- f0
  ## effect predictions 
  m <- 100
  pred_df0 <- data.frame(
    Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
    DOY = seq(min(sub_dat$DOY), max(sub_dat$DOY), length = m),
    fSite = unique(sub_dat$fSite)[1],
    Effort = 1) ## setting this to 1 net-night
  ## predictions for non random effects
  pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
  zs <- colnames(pred0$fit)
  ## for Feeagh only
  ##zs <- zs[zs != "s(trap_number):surveyIFI"]
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
    Effort = 1, ## setting this to one net
    fSite = unique(sub_dat$fSite))
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
  weight_effects_pred_1 <- rbind(weight_effects_pred_1, all_effects)
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                         DOY = 200,
                         fSite = unique(sub_dat$fSite)[1],
                         Effort = 1) ## setting this to one net
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   ##exclude = c("s(fSite)", "s(DOY)"))
                   exclude = c("s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$yhat <- exp(pred2$fit)
  pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
  pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
  weight_year_pred_1 <- rbind(weight_year_pred_1, pred_year)
}

save(weight_year_pred_1, file = "weight_year_pred.RData")
save(weight_effects_pred_1, file = "weight_effects_pred.RData")
save(weight_fits_1, file = "weight_fits.RData")
