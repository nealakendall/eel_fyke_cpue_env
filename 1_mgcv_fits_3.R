##--------------------------------
## Fit trend and covariate models
## NK: Feb 27, 2026, updated from CM
##--------------------------------

#### instead of s(year), this file uses factor(year), for the CPUE standardization ######
#### not conditioned on temperature (temp is considered a catchability, not density variable) #######

########  _yrftempc means yr as factor, watertemp as capacity covariate  ###################

### this file also contains the code to create year effects without env conditions for comparison (_noenv), with year as factor ########

##--------------------------------
## Fit trend and covariate models
##--------------------------------
library(mgcv)
library(ggplot2); theme_set(theme_bw())

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

load("all_lcdat.RData")

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

## models WITHOUT environmental characteristics, year as factor

##-----------------
## COUNT DATA FITS 
##-----------------
## containers
count_fits_noenv <- list()
effects_pred_noenv <- NULL
year_pred_noenv <- NULL


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
                         factor(Year) +
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
                         factor(Year) +
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
                         factor(Year) +
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
  gam.vcomp(f0)
  #(AIC(f0))
  #print(summary(f0))
  ##
  count_fits_noenv[[lake]] <- f0
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  # Compute means for all continuous covariates
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
      trap_depth = mean(trap_depth, na.rm = TRUE),
      trap_gradient = mean(trap_gradient, na.rm = TRUE),
    )
  ##
  years_in_model <- sort(unique(sub_dat$Year))
  keep_rows <- ! years_in_model %in% c(2025)
  ##
  pred_df2 <- data.frame(Year = factor(years_in_model[keep_rows], levels = years_in_model[keep_rows]),
                         #Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
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
  pred_year$fit_link <- pred2$fit
  pred_year$se_link  <- pred2$se.fit
  pred_year$lower_link <- pred_year$fit_link - 1.96 * pred_year$se_link
  pred_year$upper_link <- pred_year$fit_link + 1.96 * pred_year$se_link
  pred_year$yhat   <- exp(pred_year$fit_link)
  pred_year$lwr <- exp(pred_year$lower_link)
  pred_year$upr <- exp(pred_year$upper_link)
  year_pred_noenv <- rbind(year_pred_noenv, pred_year)
}

save(year_pred_noenv, file = "year_pred_noenv.RData")
save(effects_pred_noenv, file = "effects_pred_noenv.RData")
save(count_fits_noenv, file = "count_fits_noenv.RData")


#### PLOT estimated factors by year
# library(ggplot2)
# 
# ggplot(year_pred_noenv, aes(x = as.numeric(as.character(Year)), y = fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
#   geom_line(linewidth = 1.1) +
#   labs(
#     x = "Year",
#     y = "Predicted yellow eel count",
#     #title = "Yearly predicted counts with 95% CI"
#   ) +
#   theme_bw(base_size = 14) +
#   facet_wrap(~ Lake, scales = "free_y")




##-----------------
## WEIGHT DATA FITS 
##-----------------
#containers
weight_fits_noenv <- list()
weight_effects_pred_noenv <- NULL
weight_year_pred_noenv <- NULL

for(lake in lakes){
  print(lake)
  sub_dat <- subset(weightenvdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  #use this with year as a factor
  sub_dat$Year <- factor(sub_dat$Year)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake %in% c("BOH", "Bunaveela")){
    form <- as.formula(wt ~
                         factor(Year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         #guard + 
                         offset(log(Effort)))
  }
  if(lake == "Furnace"){
    form <- as.formula(wt ~
                         factor(Year) +
                         s(DOY, m = 2, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") +
                         #guard + 
                         offset(log(Effort)))
  }
  if(lake == "Feeagh"){
    form <- as.formula(wt ~
                         factor(Year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") +
                         #guard + 
                         offset(log(Effort)))
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  gam.vcomp(f0)
  #print(AIC(f0))
  #print(summary(f0))
  ###
  weight_fits_noenv[[lake]] <- f0
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  # Compute means for all continuous covariates
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
    )
  ##
  years_in_model <- sort(unique(sub_dat$Year))
  keep_rows <- ! years_in_model %in% c(2025)
  ##
  pred_df2 <- data.frame(Year = factor(years_in_model[keep_rows], levels = years_in_model[keep_rows]),
                         DOY = means$DOY,
                         fSite = unique(sub_dat$fSite)[1],
                         Effort = 1)
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   ##exclude = c("s(fSite)", "s(DOY)"))
                   exclude = c("s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$fit_link <- pred2$fit
  pred_year$se_link  <- pred2$se.fit
  pred_year$lower_link <- pred_year$fit_link - 1.96 * pred_year$se_link
  pred_year$upper_link <- pred_year$fit_link + 1.96 * pred_year$se_link
  pred_year$yhat   <- exp(pred_year$fit_link)
  pred_year$lwr <- exp(pred_year$lower_link)
  pred_year$upr <- exp(pred_year$upper_link)
  weight_year_pred_noenv <- rbind(weight_year_pred_noenv, pred_year)
}

save(weight_year_pred_noenv, file = "weight_year_pred_noenv.RData")
save(weight_effects_pred_noenv, file = "weight_effects_pred_noenv.RData")
save(weight_fits_noenv, file = "weight_fits_noenv.RData")


#### PLOT estimated factors by year
# library(ggplot2)
# 
# ggplot(weight_year_pred_noenv, aes(x = as.numeric(as.character(Year)), y = fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
#   geom_line(linewidth = 1.1) +
#   labs(
#     x = "Year",
#     y = "Predicted yellow eel weight",
#     #title = "Yearly predicted weight with 95% CI"
#   ) +
#   theme_bw(base_size = 14) +
#   facet_wrap(~ Lake, scales = "free_y")









## models WITH environmental characteristics, year as factor

##-----------------
## COUNT DATA FITS 
##-----------------
## containers
count_fits_yrftempc <- list()
effects_pred_yrftempc <- NULL
year_pred_yrftempc <- NULL


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
                         factor(Year) +
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
                         factor(Year) +
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
                         factor(Year) +
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
  gam.vcomp(f0)
  #print(AIC(f0))
  #print(summary(f0))
  ##
  count_fits_yrftempc[[lake]] <- f0
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  # Compute means for all continuous covariates
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
      trap_depth = mean(trap_depth, na.rm = TRUE),
      trap_gradient = mean(trap_gradient, na.rm = TRUE),
      watertemp = mean(watertemp, na.rm = TRUE),
      pressure = mean(pressure, na.rm = TRUE),
      wind = mean(wind, na.rm = TRUE),
      clouds = mean(clouds, na.rm = TRUE),
      moonlight = mean(moonlight, na.rm = TRUE),
      waterlev = mean(waterlev, na.rm = TRUE)
    )
  ##
  years_in_model <- sort(unique(sub_dat$Year))
  keep_rows <- ! years_in_model %in% c(2025)
  ##
  pred_df2 <- data.frame(Year = factor(years_in_model[keep_rows], levels = years_in_model[keep_rows]),
                         #Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
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
                         waterlev = mean(sub_dat$waterlev, na.rm = TRUE),
                         Effort = 2) ## 2 codends is a fyke net
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   exclude = c("s(fchain)", "s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$fit_link <- pred2$fit
  pred_year$se_link  <- pred2$se.fit
  pred_year$lower_link <- pred_year$fit_link - 1.96 * pred_year$se_link
  pred_year$upper_link <- pred_year$fit_link + 1.96 * pred_year$se_link
  pred_year$yhat   <- exp(pred_year$fit_link)
  pred_year$lwr <- exp(pred_year$lower_link)
  pred_year$upr <- exp(pred_year$upper_link)
  year_pred_yrftempc <- rbind(year_pred_yrftempc, pred_year)
}

save(year_pred_yrftempc, file = "year_pred_yrftempc.RData")
save(effects_pred_yrftempc, file = "effects_pred_yrftempc.RData")
save(count_fits_yrftempc, file = "count_fits_yrftempc.RData")


#### PLOT estimated factors by year
# library(ggplot2)
# 
# ggplot(year_pred_yrftempc, aes(x = as.numeric(as.character(Year)), y = fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
#   geom_line(linewidth = 1.1) +
#   labs(
#     x = "Year",
#     y = "Predicted yellow eel count",
#     #title = "Yearly predicted counts with 95% CI"
#   ) +
#   theme_bw(base_size = 14) +
#   facet_wrap(~ Lake, scales = "free_y")




##-----------------
## WEIGHT DATA FITS 
##-----------------
#containers
weight_fits_yrftempc <- list()
weight_effects_pred_yrftempc <- NULL
weight_year_pred_yrftempc <- NULL

for(lake in lakes){
  print(lake)
  sub_dat <- subset(weightenvdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  #use this with year as a factor
  sub_dat$Year <- factor(sub_dat$Year)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  if(lake %in% c("BOH", "Bunaveela")){
    form <- as.formula(wt ~
                         factor(Year) +
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
                         offset(log(Effort)))
  }
  if(lake == "Furnace"){
    form <- as.formula(wt ~
                         factor(Year) +
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
                         offset(log(Effort)))
  }
  if(lake == "Feeagh"){
    form <- as.formula(wt ~
                         factor(Year) +
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
                         offset(log(Effort)))
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  gam.vcomp(f0)
  ###
  weight_fits_yrftempc[[lake]] <- f0
  ##------------------------
  ## GET YEARLY PREDICTIONS
  ##------------------------
  ## here we set non-year continuous covariates to their mean and exclude random effects
  # Compute means for all continuous covariates
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
      pressure = mean(pressure, na.rm = TRUE),
      wind = mean(wind, na.rm = TRUE),
      clouds = mean(clouds, na.rm = TRUE),
      moonlight = mean(moonlight, na.rm = TRUE),
      waterlev = mean(waterlev, na.rm = TRUE)
    )
  ##
  years_in_model <- sort(unique(sub_dat$Year))
  keep_rows <- ! years_in_model %in% c(2025)
  ##
  pred_df2 <- data.frame(Year = factor(years_in_model[keep_rows], levels = years_in_model[keep_rows]),
                         DOY = means$DOY,
                         fSite = unique(sub_dat$fSite)[1],
                         watertemp = mean(sub_dat$watertemp, na.rm = TRUE),
                         pressure = means$pressure,
                         wind = means$wind,
                         clouds = means$clouds,
                         moonlight = means$moonlight,
                         waterlev = means$waterlev,
                         Effort = 1)
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                   ##exclude = c("s(fSite)", "s(DOY)"))
                   exclude = c("s(fSite)"))
  ##
  pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
  pred_year$fit_link <- pred2$fit
  pred_year$se_link  <- pred2$se.fit
  pred_year$lower_link <- pred_year$fit_link - 1.96 * pred_year$se_link
  pred_year$upper_link <- pred_year$fit_link + 1.96 * pred_year$se_link
  pred_year$yhat   <- exp(pred_year$fit_link)
  pred_year$lwr <- exp(pred_year$lower_link)
  pred_year$upr <- exp(pred_year$upper_link)
  weight_year_pred_yrftempc <- rbind(weight_year_pred_yrftempc, pred_year)
}

save(weight_year_pred_yrftempc, file = "weight_year_pred_yrftempc.RData")
save(weight_effects_pred_yrftempc, file = "weight_effects_pred_yrftempc.RData")
save(weight_fits_yrftempc, file = "weight_fits_yrftempc.RData")


#### PLOT estimated factors by year
# library(ggplot2)
# 
# ggplot(weight_year_pred_yrftempc, aes(x = as.numeric(as.character(Year)), y = fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
#   geom_line(linewidth = 1.1) +
#   labs(
#     x = "Year",
#     y = "Predicted yellow eel weight",
#     #title = "Yearly predicted weight with 95% CI"
#   ) +
#   theme_bw(base_size = 14) +
#   facet_wrap(~ Lake, scales = "free_y")

