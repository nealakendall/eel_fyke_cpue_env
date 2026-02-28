##-----------------------------
## Residual diagnostics
## NK: Feb 27, 2026, updated from CM
##-----------------------------

library(DHARMa)
library(magick)
library(mgcv)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

################ DHARMA PLOTS HERE   ##################################################################

setwd("C:\\...")

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

#count
models  <- list()
dharma  <- list()

for (lake in lakes) {
  print(lake)
  ## remove October sampling, which was out of the sampling season 
  sub_dat <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
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
  
  f0 <- gam(form,
            select = TRUE,
            method = "REML",
            family = nb(),
            data = sub_dat)
  
  models[[lake]] <- f0
  dharma[[lake]] <- simulateResiduals(f0, n = 1000)
}

res <- simulateResiduals(
  fittedModel = f0,
  n = 1000   # increase if you want more precision
)

plot(dharma[["BOH"]])
plot(dharma[["Furnace"]])
plot(dharma[["Feeagh"]])
plot(dharma[["Bunaveela"]])

for (lake in lakes) {
  png(filename = paste0("DHARMa_", lake, ".png"),
    width = 8, height = 5, units = "in", res = 300)
  
  plot(dharma[[lake]], main = "")   # no per-panel titles
  
  mtext(
    paste(lake),
    outer = TRUE,
    cex = 1.2,
    line = -1
  )
  dev.off()
}

imgs <- image_read(paste0("DHARMa_", lakes, ".png"))
row1 <- image_append(imgs[1:2], stack = FALSE)
row2 <- image_append(imgs[3:4], stack = FALSE)
combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined, "DHARMa_all_lakes_count.png")





#weight
models  <- list()
dharma  <- list()

for (lake in lakes) {
  print(lake)
  sub_dat <- subset(weightenvdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  ## sum to zero contrasts - doesn't matter currently
  contrasts(sub_dat$fSite) <- contr.sum
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
                         offset(log(Effort)))
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
                         offset(log(Effort)))
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
                         offset(log(Effort)))
  }
  ##
  f0 <- gam(form,
            select = TRUE,
            ##method = "ML",
            family = tw(),
            data = sub_dat)
  
  models[[lake]] <- f0
  dharma[[lake]] <- simulateResiduals(f0, n = 1000)
}

res <- simulateResiduals(
  fittedModel = f0,
  n = 1000   # increase if you want more precision
)

plot(dharma[["BOH"]])
plot(dharma[["Furnace"]])
plot(dharma[["Feeagh"]])
plot(dharma[["Bunaveela"]])

for (lake in lakes) {
  png(filename = paste0("DHARMa_", lake, ".png"),
      width = 8, height = 5, units = "in", res = 300)
  
  plot(dharma[[lake]], main = "")   # no per-panel titles
  
  mtext(
    paste(lake),
    outer = TRUE,
    cex = 1.2,
    line = -1
  )
  
  dev.off()
}

imgs <- image_read(paste0("DHARMa_", lakes, ".png"))
row1 <- image_append(imgs[1:2], stack = FALSE)
row2 <- image_append(imgs[3:4], stack = FALSE)
combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined, "DHARMa_all_lakes_weight.png")







######## other residual plots here ##############################################################################

setwd("C:\\...")

load("count_fits.RData")
load("all_lcdat.RData")
load("weight_fits.RData")

##------------------------
## COUNT MODEL RESIDUALS - uses negative binomial distribution
##------------------------
lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

today <- format(Sys.time(), "_%d_%m_%Y")

cex_lab <- 0.9

png(paste0("Fig_S0_residual_diagnostics", today, ".png"), height = 8, width = 7, units = "in", res = 400)
par(mfrow = c(4, 3), mar = c(2, 3, 1, 1), oma = c(2, 1, 1, 1))
set.seed(101)
for(lake in lakes){
  isb <- lake == "Bunaveela"
  f0 <- count_fits[[lake]]
  ## Calculate DHARMa randomised quantile residuals
  simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
  resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
  qqnorm(resids, bty = "l", main = "", col = "slategrey"); qqline(resids)
  legend("topleft", legend = lake, bty = "n", cex = 1.1)
  if(isb){mtext(side = 1, text = "Theoretical quantile", line = 2.5, cex = cex_lab)}
  ##
  xlim <- c(-1, 1) * max(abs(resids))
  hist(resids, probability = TRUE, bty = "l", xlim = xlim, breaks = 30, border = "lightgrey", main = "")
  curve(dnorm(x), col = "red", add = TRUE, n = 1e3)
  abline(v = 0, lty = 2)
  if(isb){mtext(side = 1, text = "Quantile residual", line = 2.5, cex = cex_lab)}
  ##
  plot(predict(f0), resids, bty = "l", col = "slategrey")
  abline(h = 0, lty = 2)
  if(isb){mtext(side = 1, text = "ln(Fitted value)", line = 2.5, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Sample quantile", line = -0.5, outer = TRUE, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Probability", line = -18, outer = TRUE, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Quantile residual", line = -35, outer = TRUE, cex = cex_lab)}
}
dev.off()


## residual autocorrelation across the fyke chain
acf_df <- NULL

##png(paste0("Fig_S2_residual_autocorrelation_", today, ".png"), height = 8, width = 7, units = "in", res = 400)
## autocorrelation
set.seed(101)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
for(lake in lakes){
  isb <- lake == "Bunaveela"
  f0 <- count_fits[[lake]]
  simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
  resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
  sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
  ##
  chains <- levels(sub_dat$fchain)
  sub_dat$resid <- resids
  ar1_vec <- sapply(chains, function(z){
    tmp <- subset(sub_dat, chain == z)
    tmp <- tmp[order(tmp$trap_number),]
    acf(tmp$resid, plot = FALSE, na.action = na.omit)$acf[2,1,1]
  })
  tmp <- data.frame(Lake = lake, fchain = chains, acf = ar1_vec)
  acf_df <- rbind(acf_df, tmp)
  ## what would a random distibution look like?
  ##rand_vec <- sapply(1:length(chains), function(z){
  rand_vec <- sapply(1:1e5 , function(z){
    y <- rnorm(20)
    ## first order unbiased - true value of mean used
    ###1/19 * sum(y[1:19] * y[2:20])
    acf(y, plot = FALSE)$acf[2,1,1]
  })
  hist(ar1_vec, breaks = seq(-1, 1, by = 0.05), xlim = c(-1, 1), border = "lightgrey", col = "grey", probability = TRUE, ylim = c(0, 3.5), main = "", xlab = "", ylab = "")
  n <- 19 ## pairs
  ##curve((1 - x^2)^((n-4)/2) / beta(a = 1/2, b = (n-2)/2), col = "blue", add = TRUE)
  lines(density(rand_vec), col = "black")
  legend("topleft", legend = lake, bty = "n")
  abline(v = 0, lty = 2)
  ##lines(density(rand_vec), col = "red")
}
mtext(side = 1, line = 0.5, text = "Trap autocorrelation", outer = TRUE)
mtext(side = 2, line = 0.5, text = "Density", outer = TRUE)
##legend("topright", legend = c("Sample acf", "Random acf", "Unbiased acf"), lty = c(NA, 1, 1), col = c("grey", "black", "blue"), pch = c(15, NA, NA), bty = "n")
legend("topright", legend = c("Observed acf", "Random acf"), lty = c(NA, 1), col = c("grey", "black"), pch = c(15, NA), bty = "n")
##dev.off()

acf_df2 <- merge(unique(lcdat[, c("Year", "Date", "fSite", "fchain")]), acf_df)

## critical bands
cb <- acf_df2[, c("Lake", "fSite")]
cb$n <- 20
cb$n[cb$fSite == "IFI"] <- 10
cb$crit <- with(cb, 1.96 / (sqrt(n -1)))

png(paste0("Figure_SX_acf_time_", today, ".png"), height = 9, width = 8, units = "in", res = 400)
ggplot(acf_df2, aes(x = Date, y = acf)) +
  geom_point() +
  facet_wrap(~ Lake + fSite) +
  geom_hline(yintercept = 0, lty = 1) +
  geom_hline(data = cb, aes(yintercept = c(1, -1) * crit), lty = 2) +
  xlab("Date") +
  ylab("First-order autocorrelation of residuals along a chain")
dev.off()

## residuals vs month

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

p_tab <- NULL

lcdat$fMonth <- factor(lcdat$Month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))

png(paste0("Fig_S3_residual_vs_month_", today, ".png"), height = 8, width = 7, units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
for(lake in lakes){
  simulationOutput <- simulateResiduals(fittedModel = count_fits[[lake]], plot = F)
  res <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
  sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)    
  boxplot(res ~ sub_dat$fMonth, notch = TRUE, bty = "l")
  legend("topleft", legend = lake, bty = "n")
  abline(h = 0, lty = 1)
}
mtext(side = 1, text = "Month", line = 0, outer = TRUE)
mtext(side = 2, text = "Residual", line = 0, outer = TRUE)
dev.off()

##------------------------
## WEIGHT MODEL RESIDUALS - uses Tweedie distribution
##------------------------
load("weight_fits.RData")
load("wdat.RData")

png(paste0("Fig_S0_residual_diagnostics_weights", today, ".png"), height = 8, width = 7, units = "in", res = 400)
set.seed(101)
par(mfrow = c(4, 3), mar = c(2, 3, 1, 1), oma = c(2, 1, 1, 1))
for(lake in lakes){
  isb <- lake == "Bunaveela"
  f0 <- weight_fits[[lake]]
  ##
  simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
  resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
  qqnorm(resids, bty = "l", main = "", col = "slategrey"); qqline(resids)
  legend("topleft", legend = lake, bty = "n", cex = 1.1)
  if(isb){mtext(side = 1, text = "Theoretical quantile", line = 2.5, cex = cex_lab)}
  ##
  xlim <- c(-1, 1) * max(abs(resids))
  hist(resids, probability = TRUE, bty = "l", xlim = xlim, breaks = 30, border = "lightgrey", main = "")
  curve(dnorm(x), col = "red", add = TRUE, n = 1e3)
  abline(v = 0, lty = 2)
  if(isb){mtext(side = 1, text = "Quantile residual", line = 2.5, cex = cex_lab)}
  ##
  plot(predict(f0), resids, bty = "l", col = "slategrey")
  abline(h = 0, lty = 2)
  if(isb){mtext(side = 1, text = "ln(Fitted value)", line = 2.5, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Sample quantile", line = -0.5, outer = TRUE, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Probability", line = -18, outer = TRUE, cex = cex_lab)}
  if(isb){mtext(side = 2, text = "Quantile residual", line = -35, outer = TRUE, cex = cex_lab)}
}
dev.off()








#join environmental indicator data to raw data used in modeling
countdat <- lcdat
#countdat <- read.csv("lcdat.csv", header=TRUE)
#countdat$newdate <- as.Date(countdat$Date, format = "%m/%d/%Y")
#countdat$DOY <- yday(countdat$newdate)

weightdat <- wdat
#weightdat <- read.csv("wdat.csv", header=TRUE)
#names(weightdat)[6] <- "DOY"

envindic <- read.csv("Env indicators by date lake.csv", header=TRUE)
#envindic$Date <- ymd(paste0(envindic$Year, "-01-01")) + days(envindic$DOY - 1)

lenenvdat <- left_join(countdat, envindic, by = c("Year", "DOY", "Lake"))
weightenvdat <- left_join(weightdat, envindic, by = c("Year", "DOY", "Lake"))

#get raw residuals
for(lake in lakes){
  f0 <- count_fits[[lake]]
  raw_residuals <- residuals(f0, type = "response")
  lenenvdat$resid[lenenvdat$Lake==lake] <- raw_residuals
}

#get raw residuals
for(lake in lakes){
  f0 <- weight_fits[[lake]]
  raw_residuals <- residuals(f0, type = "response")
  weightenvdat$resid[weightenvdat$Lake==lake] <- raw_residuals
}

#plot residuals vs. environmental variables
##LENGTH DATA
#water temp
ggplot(lenenvdat, aes(x = watertemp, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#pressure
ggplot(lenenvdat, aes(x = pressure, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#wind
ggplot(lenenvdat, aes(x = wind, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#clouds
ggplot(lenenvdat, aes(x = clouds, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#clouds:lunarlight
ggplot(lenenvdat, aes(x = cloudsvlunar, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")

##WEIGHT DATA
#water temp
ggplot(weightenvdat, aes(x = watertemp, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#pressure
ggplot(weightenvdat, aes(x = pressure, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#wind
ggplot(weightenvdat, aes(x = wind, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#clouds
ggplot(weightenvdat, aes(x = clouds, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
#clouds:lunarlight
ggplot(weightenvdat, aes(x = cloudsvlunar, y = resid)) +
  geom_point() +
  facet_wrap(~ Lake) +
  theme_minimal() +
  labs(title = "Env indicators vs residuals for each lake")
