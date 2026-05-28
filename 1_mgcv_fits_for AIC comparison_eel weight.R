############### for factor(year) #############

##--------------------------------
## Fit trend and covariate models--only for AIC comparison for different models (different env covariates included)
## CM: 5/12/22
##--------------------------------
library(mgcv)
library(ggplot2); theme_set(theme_bw())

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"


##-----------------
## WEIGHT PER EEL DATA FITS 
##-----------------

lakes <- c("Feeagh", "Furnace", "Bunaveela", "BOH")

WModelInfo <- data.frame(
  Lake = character(0),  # character vector (empty strings)
  Term = character(0),
  P_value = numeric(0),       # numeric vector (all zeros)
  AIC = numeric(0), 
  DevExpl = numeric(0),
  Model = character(0),
  ModelNum = integer(0)
  )
i <- 1

##1--MODEL WITH NO ENVIRONMENTAL COVARIATES
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
  if(lake != "Furnace"){
    form <- as.formula(wt ~
                         factor(Year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") 
                         #guard + 
                         #offset(log(Effort))
    )
  }else{
    form <- as.formula(wt ~
                         factor(Year) +
                         s(DOY, m = 2, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") 
                         #guard + 
                         #offset(log(Effort))
    )
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  smooth_terms <- rownames(summary(f0)$s.table)
  smooth_pvals <- summary(f0)$s.table[, "p-value"]
  # Combine into a data frame
  tempdat <- data.frame(
    ModelNum = i,
    Lake = lake,
    Model = "no env covariates",
    Term = smooth_terms,
    P_Value = smooth_pvals,
    DevExpl = summary(f0)$dev.expl,
    AIC = AIC(f0)
  )
  WModelInfo <- rbind(WModelInfo,tempdat)
}
  i <- i + 1
  
  ##2--MODEL WITH ENVIRONMENTAL COVARIATES
  #ALL variables as fixed effects
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
                           s(waterlev, m = 2, bs = bs_other) 
                           #guard + 
                           #offset(log(Effort))
      )
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
                           s(waterlev, m = 2, bs = bs_other) 
                           #guard + 
                           #offset(log(Effort))
      )
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
                           s(waterlev, m = 2, bs = bs_other) 
                           #guard + 
                           #offset(log(Effort))
      )
    }
    ##
    f0 <- gam(form,
              select = TRUE,
              ##method = "ML",
              family = tw(),
              data = sub_dat)
    smooth_terms <- rownames(summary(f0)$s.table)
    smooth_pvals <- summary(f0)$s.table[, "p-value"]
    # Combine into a data frame
    tempdat <- data.frame(
      ModelNum = i,
      Lake = lake,
      Model = "all env covariates",
      Term = smooth_terms,
      P_Value = smooth_pvals,
      DevExpl = summary(f0)$dev.expl,
      AIC = AIC(f0)
    )
    WModelInfo <- rbind(WModelInfo,tempdat)
  }
  i <- i + 1

#print(WModelInfo)
write.csv(WModelInfo, file = "WModelInfo.csv", row.names = TRUE)









##-----------------
## YEAR AS SMOOTH
##-----------------

library(mgcv)
library(ggplot2); theme_set(theme_bw())

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"


##-----------------
## WEIGHT DATA FITS, YEAR AS SMOOTH
##-----------------

lakes <- c("Feeagh", "Furnace", "Bunaveela", "BOH")

WModelInfo <- data.frame(
  Lake = character(0),  # character vector (empty strings)
  Term = character(0),
  P_value = numeric(0),       # numeric vector (all zeros)
  AIC = numeric(0), 
  DevExpl = numeric(0),
  Model = character(0),
  ModelNum = integer(0)
)
i <- 1

##1--MODEL WITH NO ENVIRONMENTAL COVARIATES
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
  if(lake %in% c("BOH", "Bunaveela", "Furnance")){
    form <- as.formula(wt ~
                         s(Year, m = 2, bs = bs_year) +
                         s(DOY, m = 2, bs = bs_other) +
                         s(fSite, bs = "re") 
                         #guard + 
                         #offset(log(Effort))
    )
  }else{
    form <- as.formula(wt ~
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, m = 1, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") 
                         #guard + 
                         #offset(log(Effort))
    )
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  smooth_terms <- rownames(summary(f0)$s.table)
  smooth_pvals <- summary(f0)$s.table[, "p-value"]
  # Combine into a data frame
  tempdat <- data.frame(
    ModelNum = i,
    Lake = lake,
    Model = "no env covariates",
    Term = smooth_terms,
    P_Value = smooth_pvals,
    DevExpl = summary(f0)$dev.expl,
    AIC = AIC(f0)
  )
  WModelInfo <- rbind(WModelInfo,tempdat)
}
i <- i + 1

##2--MODEL WITH ENVIRONMENTAL COVARIATES
#ALL variables as fixed effects
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
                         s(waterlev, m = 2, bs = bs_other) 
                         #guard + 
                         #offset(log(Effort))
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
                         s(waterlev, m = 2, bs = bs_other) 
                         #guard + 
                         #offset(log(Effort))
    )
  }
  if(lake == "Feeagh"){
    form <- as.formula(wt ~
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, m = 1, bs = bs_other) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 1, bs = bs_other) +
                         s(pressure, m = 1, bs = bs_other) +
                         s(wind, m = 1, bs = bs_other) +
                         s(clouds, m = 1, bs = bs_other) +
                         s(moonlight, m = 1, bs = bs_other) +
                         ti(clouds, moonlight, m = 1, bs = bs_other) +
                         s(waterlev, m = 1, bs = bs_other) 
                         #guard + 
                         #offset(log(Effort))
    )
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  smooth_terms <- rownames(summary(f0)$s.table)
  smooth_pvals <- summary(f0)$s.table[, "p-value"]
  # Combine into a data frame
  tempdat <- data.frame(
    ModelNum = i,
    Lake = lake,
    Model = "all env covariates",
    Term = smooth_terms,
    P_Value = smooth_pvals,
    DevExpl = summary(f0)$dev.expl,
    AIC = AIC(f0)
  )
  WModelInfo <- rbind(WModelInfo,tempdat)
}
i <- i + 1

#print(WModelInfo)
write.csv(WModelInfo, file = "WModelInfo.csv", row.names = TRUE)

