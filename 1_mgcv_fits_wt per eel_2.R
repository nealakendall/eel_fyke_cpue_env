##--------------------------------
## Fit trend and covariate models
## NK: Feb 27, 2026, updated from CM
##--------------------------------

#### instead of s(year), this code uses factor(year), for the CPUE standardization ######

########  We are standardizing the year effect for temperature  ###################
########  _yrftempd means yr as factor, watertemp as density covariate  ###################

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

#### WITH environmental covariates, year as factor ####

##-----------------
## WEIGHT PER EEL DATA FITS 
##-----------------
#containers
weight_fits_yrftempd <- list()
weight_effects_pred_yrftempd <- NULL
weight_year_pred_yrftempd <- NULL

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
                         s(waterlev, m = 2, bs = bs_other)) 
                         #guard + 
                         #offset(log(Effort)))
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
                         s(waterlev, m = 2, bs = bs_other)) 
                         #guard + 
                         #offset(log(Effort)))
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
  ###
  weight_fits_yrftempd[[lake]] <- f0
  ##------------------------
  ## GET YEARLY PREDICTIONS--now standardized by setting temperature to year fitted values
  ##------------------------
  ## here we set non-year and temperature continuous covariates to their mean and exclude random effects
  # 1) Compute means for all continuous covariates (except temp)
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
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
  ##
  keep_rows <- ! year_temp$Year %in% c(2025)
  #keep_rows <- ! year_temp$Year %in% c(2023, 2025)  
  ##Need to get water level data for 2022 onwards so that 2023 can be estimated well!
  
  pred_df2 <- data.frame(Year = factor(year_temp$Year[keep_rows], 
                                       levels = year_temp$Year[keep_rows]),
                         DOY = means$DOY,
                         fSite = unique(sub_dat$fSite)[1],
                         watertemp = year_temp$watertemp[keep_rows],   # <-- CONDITIONED on actual temps
                         pressure = means$pressure,
                         wind = means$wind,
                         clouds = means$clouds,
                         moonlight = means$moonlight,
                         waterlev = means$waterlev)
                         #Effort = 1)
  ##
  pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
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
  weight_year_pred_yrftempd <- rbind(weight_year_pred_yrftempd, pred_year)
}

save(weight_year_pred_yrftempd, file = "weight_year_pred_yrftempd.RData")
save(weight_effects_pred_yrftempd, file = "weight_effects_pred_yrftempd.RData")
save(weight_fits_yrftempd, file = "weight_fits_yrftempd.RData")
