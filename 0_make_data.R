##--------------------------------
## Make the datasets for the models
## NK: Feb 27, 2026, updated from CM
##--------------------------------
library(gdata)
library(reshape)  #install.packages("reshape")  #if needed
library(ggplot2); theme_set(theme_bw())
library(lubridate)  #install.packages("lubridate") #if needed
library(dplyr)
library(tidyr)

setwd("C:\\....")
rm(list=ls())

##------------
## COUNT DATA
##------------
dat <- read.csv("Fyke surveys.csv", header=TRUE)

## remove October data
dat <- subset(dat, Month != "Oct")

## merge in the depths
depths <- read.csv("Sites Depths_from Coilin.csv")

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

names(depths)[names(depths) == "Inshore.Dept"] <- "dstart"
names(depths)[names(depths) == "Dept.at.75m...midway"] <- "dmid"
names(depths)[names(depths) == "Depth.at.outer.End"] <- "dfinish"

depths$Lake <- depths$Location

## don't know which way the chain oriented for IFI data - assuming average depth

## first fill in missing outer depths for IFI data
idx <- which(depths$Comment == "IFI random survey" & is.na(depths$X2009.2010.Depth.finish..m.))
depths$X2009.2010.Depth.finish..m.[idx] <- depths$X2009..2010.Depth.start..m.[idx]

## get average depth
idx <- which(depths$Comment == "IFI random survey")

depths$dstart[idx] <- depths$dmid[idx] <- depths$dfinish[idx] <-
  apply(depths[idx, c("X2009..2010.Depth.start..m.", "X2009.2010.Depth.finish..m.")], 1, function(x){sum(x)/2})

## merge in the 3 point depths
dat <- merge(dat, depths[, c("Lake", "Site", "dstart", "dmid", "dfinish")], all.x = TRUE)

idx <- which(is.na(dat$dstart) | is.na(dat$dmid) | is.na(dat$dfinish))

## remove those with no depths
## dat[idx, c("Lake", "Site")]
dat <- dat[-idx,] 

names(dat)[names(dat) == "Num.of.Nights"] <- "Effort"

## codend data
xs <- paste0("X", 1:20)
cdat <- dat[, c("Year", "Month", "Date", "Lake", "Code", "Effort", "dstart", "dmid", "dfinish", "Site", xs)]

## long form
lcdat <- melt(cdat, id.vars = c("Year", "Month", "Date", "Lake", "Code", "Effort", "dstart", "dmid", "dfinish", "Site"))

lcdat$trap_number <- as.numeric(gsub("X", "", as.character(lcdat$variable)))

## remove NA observations
lcdat <- subset(lcdat, !is.na(value))

## get trap-specific linear depths
## uses trap mouth position in metres along chain
## see: data/Trap_depths_along_chain.xlsx

get_dg <- function(tn, d0, d1, d2){
  ##------------
  ## returns inferred depth and gradient at trap mouth
  ## tn: trap number
  ## d0: starting depth (metres)
  ## d1: mid depth (metres)
  ## d2: finishing depth (metres)
  ##------------
  ## distance to trap mouth 
  tdist <- cumsum(c(3.4, rep(c(8, 7), times = 10)[-20]))
  names(tdist) <- 20:1 ## trap 1 corresponds to depth d2
  b0 <- d0
  mid <- 75
  ## slope in first section
  b1 <- (d1 - d0) / mid
  ## slope in second section
  b2 <- (d2 - d1) / mid
  x <- tdist[as.character(tn)]
  idx <- x <= mid
  depth <- gradient <- rep(NA, length(tn))
  depth[idx] <- (b0 + b1 * x)[idx]
  depth[!idx] <- (b0 + b1 * mid + b2 * (x - mid))[!idx]
  ## gradient
  gradient[idx] <- b1[idx]
  gradient[!idx] <- b2[!idx]
  res <- list(depth = depth, gradient = gradient)
  return(res)
}

dg <- with(lcdat, get_dg(tn = trap_number, d0 = dstart, d1 = dmid, d2 = dfinish))

## correct for Bunaveela 2009, 2010
## only 5 pairs fished
idx <- which(lcdat$Lake == "Bunaveela" & lcdat$Year %in% c(2009, 2010))

tdist <- cumsum(c(3.4, rep(c(8, 7), times = 5)[-10]))
names(tdist) <- 10:1 
b0 <- lcdat$dstart[idx]
mid <- 75
## slope in first section
b1 <- (lcdat$dmid[idx] - lcdat$dstart[idx]) / mid
x <- tdist[as.character(lcdat$trap_number[idx])]
dg$depth[idx] <- (b0 + b1 * x)
dg$gradient[idx] <- b1

## changing depths to negative
lcdat$trap_depth <- -dg$depth
lcdat$trap_gradient <- dg$gradient

idx <- which(lcdat$trap_number %in% seq(2, 20, by = 2))
lcdat$trap_gradient[idx] <- -lcdat$trap_gradient[idx]
## note: negative gradients face downwards

##pdf("Bunaveela_depths.pdf", height = 6, width = 8)
##ggplot(subset(lcdat, Lake == "Bunaveela"), aes(x = trap_number, y = trap_depth, colour = Site)) +
##    geom_point() + facet_wrap(~Year)
##dev.off()

# lakes <- unique(lcdat$Lake)
# today <- format(Sys.time(), "_%d_%m_%Y")
# pdf(paste0("trap_depths", today, ".pdf"), height = 7, width = 10)
# for(l in lakes){
#   tmp <- subset(lcdat, Lake == l)
#   if(l == "Feeagh"){
#     ## remove IFI sites for depth plots
#     idx <- grep("Fee", tmp$Site)
#     tmp <- tmp[-idx,]
#     tmp <- droplevels(tmp)
#   }
#   p <- ggplot(tmp, aes(x = trap_number, y = trap_depth)) +
#     geom_point(aes(colour = trap_gradient), size = 2) +
#     facet_wrap( ~ Site) +
#     paletteer::scale_colour_paletteer_c(palette = "viridis::plasma", name = "Trap gradient") +
#     ##scale_colour_gradient(low = "yellow", high = "darkred") +
#     xlab("Trap number") +
#     ylab("Trap opening depth (m)") +
#     ggtitle(l)
#   print(p)
# }
# dev.off()

lcdat$chain <- paste0("chain", lcdat$Code)

pairs <- data.frame(X = paste0("X", 1:20), pair = rep(paste0("pair", 1:10), each = 2))

lcdat$pair <-
  pairs$pair[match(as.character(lcdat$variable), pairs$X)]

lcdat <- droplevels(lcdat)

## NB ORDER SITES according to ABC
lcdat$Site[grep("Fee", lcdat$Site)] <- "IFI"
uniq_sites <- unique(lcdat$Site)
site_levels <- toupper(letters[1:10])
site_levels <- c(site_levels, uniq_sites[!uniq_sites %in% site_levels])

lcdat$fSite <- factor(lcdat$Site, levels = site_levels)
lcdat$fchain <- factor(lcdat$chain) 
lcdat$fpair <- factor(lcdat$pair)

lcdat$fchain_fpair <- factor(with(lcdat, paste(fchain, fpair, sep = "-")))

## get the julian day
lcdat$DateNew <- trimws(lcdat$Date)
lcdat$DateNew <- as.character(lcdat$DateNew)
lcdat$DateNew <- mdy(lcdat$DateNew)
lcdat$DOY <- as.numeric(format(lcdat$DateNew, "%j"))

names(lcdat)[names(lcdat) == "value"] <- "count"

lcdat$survey <- factor(ifelse(lcdat$Site == "IFI", "IFI", "Russell"))

## trap level ID
lcdat$fID <- factor(apply(lcdat[, c("chain", "trap_number")], 1, paste, collapse = ":"))

## otter guard on/off
lcdat$guard <- ifelse(lcdat$Year >= 2015, "yes", "no")

save(lcdat, file = "all_lcdat.RData")

#write.csv(lcdat, file = "lcdat.csv", row.names = TRUE)

##-------------
## WEIGHT DATA 
##-------------
vars2keep <- c("Lake", "Site", "Year", "Month", "Julian.Day", "Total.Weight..kg.", "Effort..net.nights.", "dmid")

wdat <- dat[, vars2keep]
names(wdat)[names(wdat) == "Julian.Day"] <- "DOY"
names(wdat)[names(wdat) == "Total.Weight..kg."] <- "wt"
names(wdat)[names(wdat) == "Effort..net.nights."] <- "Effort"
names(wdat)[names(wdat) == "dmid"] <- "depth_mid"

ggplot(wdat, aes(x = Year, y = wt)) +
  geom_point() +
  facet_wrap(~ Lake, scales = "free_y")

wdat <- subset(wdat, !is.na(wt))

wdat$Site[grep("Fee", wdat$Site)] <- "IFI"
uniq_sites <- unique(wdat$Site)
site_levels <- toupper(letters[1:10])
site_levels <- c(site_levels, uniq_sites[!uniq_sites %in% site_levels])

wdat$fSite <- factor(wdat$Site, levels = site_levels)

## otter guard on/off
wdat$guard <- ifelse(wdat$Year >= 2015, "yes", "no")

save(wdat, file = "wdat.RData")

#write.csv(wdat, file = "wdat.csv", row.names = TRUE)





#####################################################################
##### Datasets with environmental indicators ###############

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

countenvdat <- left_join(countdat, envindic, by = c("Year", "DOY", "Lake"))
weightenvdat <- left_join(weightdat, envindic, by = c("Year", "DOY", "Lake"))

#write.csv(countenvdat, file = "Count env data.csv", row.names = TRUE)


### Standardize the continuous variables, by lake ###

lakes <- c("Feeagh", "Furnace", "Bunaveela", "BOH")
vars <- c("trap_depth", "trap_gradient", "watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")

for(lake in lakes){
  #print(lake)
    for(var in vars){
      #print(var)
      mean_use <- mean(countenvdat[[var]][countenvdat$Lake==lake], na.rm = TRUE)
      stdev_use <- sd(countenvdat[[var]][countenvdat$Lake==lake], na.rm = TRUE)
      rows_to_standardize <- countenvdat$Lake==lake
      x <- countenvdat[[var]][rows_to_standardize]
      standardized <- (x - mean_use) / stdev_use
      countenvdat[[var]][rows_to_standardize] <- standardized
    }   
}

vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")

for(lake in lakes){
  #print(lake)
  for(var in vars){
    #print(var)
    mean_use <- mean(weightenvdat[[var]][weightenvdat$Lake==lake], na.rm = TRUE)
    stdev_use <- sd(weightenvdat[[var]][weightenvdat$Lake==lake], na.rm = TRUE)
    rows_to_standardize <- weightenvdat$Lake==lake
    x <- weightenvdat[[var]][rows_to_standardize]
    standardized <- (x - mean_use) / stdev_use
    weightenvdat[[var]][rows_to_standardize] <- standardized
  }   
}


################################
## plot the not-standardized env covariates by lake and doy (across all years)

#count data

dat_long <- countenvdat %>%
  pivot_longer(
    cols = c(watertemp, pressure, wind, waterlev, clouds, moonlight, dmid, count),
    names_to = "covariate",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  mutate(
    covariate = factor(
      covariate,
      levels = c("watertemp", "pressure", "wind", "waterlev", "clouds", "moonlight", "dmid", "count"),
      labels = c("Water temp", "Pressure", "Wind speed", "Water level", "Cloud cover", "Moonlight", "Trap depth", "Eel count")
    ),
    Lake = factor(
      Lake,
      levels = c("BOH", "Furnace", "Feeagh", "Bunaveela")  # replace with your lake names in desired order
    )
  )

p <- ggplot(
  dat_long,
  aes(x = DOY, y = value)
) +
  geom_point(alpha = 0.4, size = 1) +
  facet_grid(covariate ~ Lake, scales = "free_y") +
  labs(
    x = "Day of year",
    y = "Environmental variable"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

# ggsave("Env conditions in surveys by doy_count.png", plot = p, width = 8, height = 10, dpi = 300)



#weight data

dat_long <- weightenvdat %>%
  pivot_longer(
    cols = c(watertemp, pressure, wind, waterlev, clouds, moonlight, wt),
    names_to = "covariate",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  mutate(
    covariate = factor(
      covariate,
      levels = c("watertemp", "pressure", "wind", "waterlev", "clouds", "moonlight", "wt"),
      labels = c("Water temp", "Pressure", "Wind speed", "Water level", "Cloud cover", "Moonlight", "Eel weight"),
    ),
    Lake = factor(
      Lake,
      levels = c("BOH", "Furnace", "Feeagh", "Bunaveela")  # replace with your lake names in desired order
    )
  )

p <- ggplot(
  dat_long,
  aes(x = DOY, y = value)
) +
  geom_point(alpha = 0.4, size = 1) +
  facet_grid(covariate ~ Lake, scales = "free_y") +
  labs(
    x = "Day of year",
    y = "Environmental variable"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

# ggsave("Env conditions in surveys by doy_weight.png", plot = p, width = 8, height = 10, dpi = 300)
