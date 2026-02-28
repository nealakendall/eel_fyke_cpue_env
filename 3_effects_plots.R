##--------------------------------
## Plots for the paper
## NK: Feb 27, 2026, updated from CM
##--------------------------------

#####  in the effects plots, use smoothed year (as I can't figure out how to do with year as factor), use _yrstempc ####
#######  in the effects plots, removed year, so now only plotted year via year effect  ####
#### in the year effects plots (CPUE standardization), instead of s(year), use factor(year)  ######

##--------------------------------
## Plots for the paper
##--------------------------------

library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(dplyr)
library(MASS) ## for random multivariate normal simulation
library(forcats)
library(ggh4x)

setwd("C:\\Users\\kendanwk7\\OneDrive - Washington State Executive Branch Agencies\\Desktop\\Burrishoole eel abundance\\Eel model_std_cloundslightinteract")

today <- format(Sys.time(), "_%d_%m_%Y")

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

#load("all_lcdat.RData")

## basic data plots
p0 <- ggplot(countenvdat, aes(x = jitter(Year), y = count)) +
  geom_point() +
  facet_wrap(~ Lake, ncol = 1) +
  ylab("Count per trap") +
  theme(strip.background = element_rect(fill = "grey"))

pdf("raw_data.pdf", height = 8, width = 10)
print(p0)
dev.off()


site_df <- as.data.frame(with(subset(countenvdat, survey != "IFI"), table(Year, Site, Lake)))

site_df <- subset(site_df, Freq > 0)

site_df$Year <- as.numeric(as.character(site_df$Year))

library(viridis)
p1 <- ggplot(site_df, aes(x = Year, y = Site, fill = Freq)) +
  geom_tile(color="white", linewidth=0.1) +   #changed "size = 0.1" to "linewidth = 0.1"
  ##facet_wrap(~ Lake, scales = "free")
  facet_wrap(~ Lake, scales="free_y", ncol = 1) +
  scale_fill_viridis(name="# Codends") +
  theme(strip.background = element_rect(fill = "grey"),
        legend.position = "none")

p2 <- ggplot(countenvdat, aes(x = Year, y = trap_depth)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~ Lake, ncol = 1) +
  ylab("Trap depth") +
  theme(strip.background = element_rect(fill = "grey"))

p3 <- ggplot(countenvdat, aes(x = Year, y = trap_gradient)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~ Lake, ncol = 1) +
  ylab("Trap gradient") +
  theme(strip.background = element_rect(fill = "grey"))

png(paste0("Variable_plots", today, ".png"), height = 6, width = 10, units = "in", res = 400)
grid.arrange(p0, p1, p2, p3, nrow = 1)
dev.off()

##---------------
## RESULTS PLOTS
##---------------

##---------------
## COUNT PLOTS
##---------------

##-----------------------
##  COVARIATE EFFECTS PLOTS
##-----------------------
load("effects_pred_yrstempc.RData")
load("count_fits_yrstempc.RData")

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

effects_pred_yrstempc$Lake <- factor(effects_pred_yrstempc$Lake, levels = lakes)

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(count_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

vars <- unique(effects_pred_yrstempc$variable)

vars <- vars[!vars %in% c("Year", "fSite", "fchain")]

for(v in vars){
  sub <- subset(effects_pred_yrstempc, variable == v)
  sub$x <- as.numeric(sub$x)
  ##sub$Lake <- factor(sub$Lake, levels = c("BOH", "Furnace", "Feeagh", "Bunaveela"))
  var_lab <- stringr::str_to_title(gsub("\\_", " ", v))
  ##
  if(v == "trap_number"){
      pv <- pvalue_df[grep(v, pvalue_df$var),]
  } else {
      smooth_name <- paste0("s(", v, ")")
      pv <- subset(pvalue_df, var == smooth_name)
      }
  pv$survey <- "Russell"
  pv$survey[grep("surveyIFI", pv$var)] <- "IFI"
  ##if(length() > 0){
  ##    pv <- pv[-grep("surveyIFI", pv$var), ] ## double survey trap number effect for IFI 
  ##}
  if(v == "trap_number"){
    tmpF <- subset(pv, Lake == "Feeagh")
    p_string <- paste(tmpF$p_value2, collapse = ", ") ## IFI first
    tmp <- subset(pv, Lake != "Feeagh")
    tmp0 <- subset(tmpF, var == "s(trap_number):surveyRussell")
    tmp0$p_value2 <- p_string
    pv <- rbind(tmp, tmp0)
  }
  p <-
    ggplot(sub, aes(x = x, y = yhat, group = survey)) +
    geom_ribbon(aes(ymin = ylwr, ymax = yupr), fill = "grey", alpha = 0.7) +
    geom_line(colour = "black", lwd = 0.5) +
    #geom_rug(sides = "b", alpha = 1/2) +  # Add rug marks at the bottom
    geom_rug(sides = "b", alpha = 1/2, position = "jitter") +  # Add rug marks at the bottom
    facet_wrap(~Lake, nrow = 1) +
    ##xlab(var_lab) +
    xlab("") +
    ylab(paste(var_lab, "effect")) +
    geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
              hjust = 1, vjust  = 1.5, size = 3) +
    theme(plot.margin = unit(c(0, 1, 0, 1), "lines"), legend.position = "none") +
    scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))
  if(v != "watertemp"){
    p <- p + theme(##strip.background = element_rect("grey"),
      strip.text.x = element_blank() , 
      strip.background = element_blank())
  }
  assign(paste0("p_", v), p)
}

## chain effects
sub <- subset(effects_pred_yrstempc, variable == "fchain")
var_lab <- "Chain effect"
pv <- pvalue_df[grep("chain", pvalue_df$var),]

p_chain <-
  ggplot(sub, aes(x = exp(yhat))) +
  geom_histogram(bins = 30, fill = "slategrey", colour = "grey") +
  facet_wrap(~Lake, nrow = 1) +
  xlab("") +
  ylab("Chain effect") +
  theme(plot.margin = unit(c(0, 1, 0, 1), "lines"),
        strip.text.x = element_blank() , 
        strip.background = element_blank(),
        legend.position = "none")  +
  geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
            hjust = 1.2, vjust = 2, size = 3) +
  scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))

## site effects
sub <- subset(effects_pred_yrstempc, variable == "fSite")
sub$x[sub$x == "Back Weir"] <- "BW"
sub$x[sub$x == "Front Weir"] <- "FW"
sub$x[sub$x == "North"] <- "N"
sub$x[sub$x == "South"] <- "S"
sub$x[sub$x == "Weir"] <- "W"

var_lab <- "Site effect"

pv <- pvalue_df[grep("Site", pvalue_df$var),]

p_site <-
  ggplot(sub, aes(y = yhat, x = x)) +
  geom_point() +
  geom_segment(aes(y = ylwr, yend = yupr, xend = x)) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~Lake, nrow = 1, scales = "free_x") +
  xlab("") +
  ylab("Site effect") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0, 1, 0, 1), "lines"),
        strip.text.x = element_blank() , 
        strip.background = element_blank(),
        legend.position = "none"          
  )  +
  geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
            hjust = 1.3, vjust = 2, size = 3) +
  scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))

##png(paste0("Count_effect_plots", today, ".png"), height = 10, width = 7, units = "in", res = 400)
jpeg(paste0("Fig3_count_effect_plots", today, ".jpg"), height = 8.5, width = 7, units = "in", res = 600)
#jpeg(paste0("FigS5_count_effect_plots_varsnotstd", today, ".jpg"), height = 8.5, width = 7, units = "in", res = 600)
#grid.arrange(p_Year, p_DOY, p_trap_depth, p_trap_gradient, p_trap_number, p_watertemp, p_pressure, p_wind, p_clouds, p_chain, p_site, ncol = 1)
grid.arrange(p_watertemp, p_pressure, p_wind, p_waterlev, p_clouds, p_moonlight, ncol = 1)  #p_waterlev having issues with the jittering
dev.off()

##-----------------------
##  YEARLY TREND PLOTS
##-----------------------
load("year_pred_noenv.RData")  #this is with year as factor and no environmental data
load("year_pred_yrftempd.RData")  #this is with year as factor and year conditioned on watertemp (watertemp as density)
load("year_pred_yrftempc.RData")  #this is with year as factor and year not conditioned on watertemp (watertemp as catchability)
load("count_fits_yrftempd.RData")

year_pred_noenv$Lake <- factor(year_pred_noenv$Lake, levels = lakes)
var_lab <- "Year effect"

countenvdat$Lake <- factor(countenvdat$Lake, levels = lakes)

mean_count <- aggregate(count ~ Year + Lake, mean, data = countenvdat)

## get the count per pair
pair_tab <- table(countenvdat$fchain_fpair)
pair_remove <- names(pair_tab)[pair_tab < 2]
## 2 unpaired codends out of 5302 pairs

pair_count <- aggregate(count ~ Year + Lake + fchain_fpair, sum, data = countenvdat)

pair_count$Lake <- fct_relevel(pair_count$Lake,
                               "BOH", "Furnace", "Feeagh", "Bunaveela")
year_pred_noenv$Lake    <- fct_relevel(year_pred_noenv$Lake,    levels(pair_count$Lake))
year_pred_yrftempd$Lake <- fct_relevel(year_pred_yrftempd$Lake, levels(pair_count$Lake))
year_pred_yrftempc$Lake <- fct_relevel(year_pred_yrftempc$Lake, levels(pair_count$Lake))

# write.csv(year_pred_noenv, file = "year_pred_noenv.csv", row.names = FALSE)
# write.csv(year_pred_yrftempd, file = "year_pred_yrftempd.csv", row.names = FALSE)
# write.csv(year_pred_yrftempc, file = "year_pred_yrftempc.csv", row.names = FALSE)

all_years <- 1987:2025

countenvdat$Year_num        <- as.numeric(as.character(countenvdat$Year))
pair_count$Year_num         <- as.numeric(as.character(pair_count$Year))
year_pred_noenv$Year_num    <- as.numeric(as.character(year_pred_noenv$Year))
year_pred_yrftempd$Year_num <- as.numeric(as.character(year_pred_yrftempd$Year))
year_pred_yrftempc$Year_num <- as.numeric(as.character(year_pred_yrftempc$Year))

off_box  <- -0.30   # boxplot + red mean CI
off_grn  <- 0       # green diamond
off_blue <- +0.30   # blue triangles

py0 <-
  ggplot()+
  facet_wrap(~Lake, scales = "free", ncol = 2) +
  # first the raw (nominal) data
  #geom_boxplot(data = pair_count, 
  #             aes(x = Year_num + off_box, y = count, group = Year)) +
  stat_summary(data = pair_count,
               aes(x = Year_num + off_box, y = count),
               fun.data = "mean_cl_boot",
               colour = "red",
               size = 0.3) +
  # then the values without environmental variables (green star)
  geom_errorbar(data = year_pred_noenv,
                aes(x = Year_num + off_grn, ymin = lwr, ymax = upr, y = yhat),
                width = 0.1, colour = "green") +
  geom_point(data = year_pred_noenv,
             aes(x = Year_num + off_grn, y = yhat),
             size = 2, shape = 18, colour = "green") +  # green diamond
  # then the values with env variables, conditioned on temp (blue triangles)
  geom_errorbar(data = year_pred_yrftempd,
                aes(x = Year_num + off_blue, ymin = lwr, ymax = upr, y = yhat),
                width = 0.1, colour = "blue") +
  geom_point(data = year_pred_yrftempd,
             aes(x = Year_num + off_blue, y = yhat),
             size = 1.5, shape = 17, colour = "blue") +  # blue triangle
  # then the values with env variables, NOT conditioned on temp (red squares)
  # geom_errorbar(data = year_pred_yrftempc,
  #               aes(x = Year_num + off_pur, ymin = lwr, ymax = upr, y = yhat),
  #               width = 0.12, colour = "purple") +
  # geom_point(data = year_pred_yrftempc,
  #            aes(x = Year_num + off_pur, y = yhat),
  #            size = 2, shape = 15, colour = "purple") +   #purple square
  ## Labels + theme
  # scale_x_continuous(
  #   breaks = seq(1987, 2025, 5),
  #   labels = seq(1987, 2025, 5),
  #   limits = c(min(all_years)-0.5, max(all_years)+0.5)) +
  xlab("Year") +
  ylab("Count per fyke net (2 codends)") +
  theme(strip.background = element_rect(fill = "grey"))

jpeg(paste0("Fig5_count_per_year_plots", today, ".jpg"), height = 4.5, width = 10, units = "in", res = 600)
py0
dev.off()


#PLOT comparing raw data to standardized value (conditioned on temp)
y_scales <- list(
  BOH        = scale_y_continuous(limits = c(0, 80)),
  Furnace   = scale_y_continuous(limits = c(0, 20)),
  Feeagh    = scale_y_continuous(limits = c(0, 8)),
  Bunaveela = scale_y_continuous(limits = c(0, 6))
)

#get the right order of the lakes
lake_levels <- c("BOH", "Furnace", "Feeagh", "Bunaveela")
pair_count2 <- pair_count %>%
  mutate(Lake = factor(Lake, levels = lake_levels))
year_pred2 <- year_pred_yrftempd %>%
  mutate(Lake = factor(Lake, levels = lake_levels))

py1 <-
  ggplot() +
  facet_wrap(~Lake, ncol = 4, scales = "free_y") +
  geom_boxplot(
    data = pair_count2,
    aes(x = Year_num, y = count, group = Year)
  ) +
  stat_summary(
    data = pair_count2,
    aes(x = Year_num, y = count),
    fun.data = "mean_cl_boot",
    colour = "red",
    size = 0.3
  ) +
  geom_ribbon(
    data = year_pred2,
    aes(x = Year_num, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.40,
    fill = "lightblue"
  ) +
  geom_line(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    colour = "blue",
    linewidth = 0.5
  ) +
  geom_point(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    size = 2, shape = 15, colour = "blue"
  ) +
  facetted_pos_scales(
    y = y_scales
  ) +
  xlab("Year") +
  ylab("Count per fyke net (2 codends)") +
  theme(strip.background = element_rect(fill = "grey"))


py2 <-
  ggplot() +
  facet_wrap(~Lake, ncol = 4, scales = "free_y") +
  geom_ribbon(
    data = year_pred2,
    aes(x = Year_num, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.40,
    fill = "lightblue"
  ) +
  geom_line(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    colour = "blue",
    linewidth = 0.5
  ) +
  geom_point(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    size = 2, shape = 15, colour = "blue"
  ) +
  xlab("Year") +
  ylab("Count per fyke net (2 codends)") +
  theme(strip.background = element_rect(fill = "grey"))



## estimate and plot percentage change between the start and the end of the time series
percent_df <- NULL

for(lake in lakes){
  f0 <- count_fits_yrftempd[[lake]]
  sub_dat <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
  # 1) Compute means for all continuous covariates (except temp)
  means <- sub_dat |> 
    summarise(
      DOY = mean(DOY, na.rm = TRUE),
      trap_depth = mean(trap_depth, na.rm = TRUE),
      trap_gradient = mean(trap_gradient, na.rm = TRUE),
      pressure = mean(pressure, na.rm = TRUE),
      wind = mean(wind, na.rm = TRUE),
      clouds = mean(clouds, na.rm = TRUE),
      moonlight = mean(moonlight, na.rm = TRUE),
      waterlev = mean(waterlev, na.rm = TRUE)
    )
  # 2) Compute YEAR-SPECIFIC mean temperature
  year_temp <- sub_dat |> 
    group_by(Year) |> 
    summarise(watertemp = mean(watertemp, na.rm = TRUE))
  ##
  pred_df2 <- data.frame(Year = year_temp$Year,
                         DOY = means$DOY,
                         trap_depth = means$trap_depth, 
                         trap_gradient = means$trap_gradient,
                         pressure = means$pressure,
                         wind = means$wind,
                         clouds = means$clouds,
                         moonlight = means$moonlight,
                         waterlev = means$waterlev,
                         watertemp = year_temp$watertemp,  # <-- CONDITIONED on actual temps
                         fSite = unique(sub_dat$fSite)[1],
                         fchain = unique(sub_dat$fchain)[1],
                         trap_number = 10,
                         survey = "Russell",
                         Effort = 2)
  ##
  if(lake == "Feeagh"){
    Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fchain)", "s(fSite)"))
  }else{
    Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fchain)", "s(fSite)"))
  }
  ## draw from posterior - see: https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
  br <- mvrnorm(10000,coef(f0), vcov(f0))
  lnN <- t(Xp %*% t(br))
  ## percentage decline
  #percent_decline <- 100 * (1 - exp(lnN[, 2] - lnN[, 1]))
  if(lake == "BOH"){
          percent_decline <- 100 * (1 - exp(lnN[, 15] - lnN[, 1]))  #because data for 2025 not available
          }
  if(lake == "Bunaveela"){
    percent_decline <- 100 * (1 - exp(lnN[, 21] - lnN[, 1])) #because data for 2025 not available
          }
  if(lake == "Furnace"){
    percent_decline <- 100 * (1 - exp(lnN[, 19] - lnN[, 1])) #because data for 2025 not available
          }
  if(lake == "Feeagh"){
    percent_decline <- 100 * (1 - exp(lnN[, 22] - lnN[, 1])) #because data for 2025 not available
          }
  df <- data.frame(Lake = lake, pd = percent_decline)
  percent_df <- rbind(percent_df, df)
  rm(df)
}

percent_df$Lake <- factor(percent_df$Lake, levels = c("BOH", "Furnace", "Feeagh", "Bunaveela"))

agg0 <- aggregate(pd ~ Lake, FUN = mean, data = percent_df)
agg1 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.025, data = percent_df)
names(agg1)[2] <- c("lwr")
agg2 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.975, data = percent_df)
names(agg2)[2] <- c("upr")

summary_df <- merge(merge(agg0, agg1), agg2)

write.csv(summary_df, file = "count_percent_decline.csv", row.names = FALSE)

#specify lake order for the plots:
lake_order <- c("BOH", "Furnace", "Feeagh", "Bunaveela")
percent_df$Lake <- factor(percent_df$Lake, levels = lake_order)
summary_df$Lake <- factor(summary_df$Lake, levels = lake_order)

py3 <- ggplot(percent_df, aes(x = pd)) +
  geom_density(fill = "grey") +
  facet_wrap(~ Lake, ncol = 4, scales = "free_y") +
  ylab("Density") +
  geom_vline(data = summary_df, aes(xintercept = pd)) +
  geom_vline(data = summary_df, aes(xintercept = lwr), lty = 3) +
  geom_vline(data = summary_df, aes(xintercept = upr), lty = 3) +
  geom_vline(xintercept = 0, lty = 4) +
  xlab("Percentage decline in standardised count trend 1987-2025") +
  coord_cartesian(xlim = c(-10, 100))

##png(paste0("Count_percent_decline_plots", today, ".png"), height = 8, width = 10, units = "in", res = 400)
jpeg(paste0("Fig6_count_percent_decline_plots", today, ".jpg"), height = 7, width = 10, units = "in", res = 600)
grid.arrange(py1, py2, py3, nrow =3)
dev.off()




##---------------
## WEIGHTS PLOTS
##---------------

##-----------------------
##  COVARIATE EFFECTS PLOTS
##-----------------------
load("weight_effects_pred_yrstempc.RData")
load("weight_fits_yrstempc.RData")
#load("wdat.RData")

weightenvdat <- subset(weightenvdat, Month != "Oct" & !is.na(wt))

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

weight_effects_pred_yrstempc$Lake <- factor(weight_effects_pred_yrstempc$Lake, levels = lakes)

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(weight_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

vars <- unique(weight_effects_pred_yrstempc$variable)

vars <- vars[!vars %in% c("fSite", "fchain")]

for(v in vars){
  sub <- subset(weight_effects_pred_yrstempc, variable == v)
  sub$x <- as.numeric(sub$x)
  var_lab <- stringr::str_to_title(gsub("\\_", " ", v))
  ##
  smooth_name <- paste0("s(", v, ")")
  pv <- subset(pvalue_df, var == smooth_name)
  #pv <- pvalue_df[grep(v, pvalue_df$var),]
  if(length(grep("surveyIFI", pv$var)) > 0){
    pv <- pv[-grep("surveyIFI", pv$var), ] ## double survey trap number effect for IFI 
  }
  p <- ggplot(sub, aes(x = x, y = yhat)) +
    geom_ribbon(aes(ymin = ylwr, ymax = yupr), fill = "grey", alpha = 0.7) +
    geom_line(colour = "black", lwd = 0.5) +
    #geom_rug(sides = "b", alpha = 1/2) +  # Add rug marks at the bottom
    geom_rug(sides = "b", alpha = 1/2, position = "jitter") +  # Add rug marks at the bottom
    facet_wrap(~Lake, nrow = 1) +
    ##xlab(var_lab) +
    xlab("") +
    ylab(paste(var_lab, "effect")) +
    geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
              hjust = 1, vjust  = 1.5, size = 3) +
    theme(plot.margin = unit(c(0, 1, 0, 1), "lines"), legend.position = "none") +
    scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))
  if(v != "watertemp"){
    p <- p + theme(##strip.background = element_rect("grey"),
      strip.text.x = element_blank() , 
      strip.background = element_blank())
  }
  assign(paste0("p_", v), p)
}

## site effects
sub <- subset(weight_effects_pred_yrstempc, variable == "fSite")
sub$x[sub$x == "Back Weir"] <- "BW"
sub$x[sub$x == "Front Weir"] <- "FW"
sub$x[sub$x == "North"] <- "N"
sub$x[sub$x == "South"] <- "S"
sub$x[sub$x == "Weir"] <- "W"

var_lab <- "Site effect"

pv <- pvalue_df[grep("Site", pvalue_df$var),]

p_site <-
  ggplot(sub, aes(y = yhat, x = x)) +
  geom_point() +
  geom_segment(aes(y = ylwr, yend = yupr, xend = x)) +
  geom_hline(yintercept = 0, lty = 2) +
  facet_wrap(~Lake, nrow = 1, scales = "free_x") +
  xlab("") +
  ylab("Site effect") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(0, 1, 0, 1), "lines"),
        strip.text.x = element_blank() , 
        strip.background = element_blank(),
        legend.position = "none"          
  )  +
  geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
            hjust = 1.3, vjust = 2, size = 3) +
  scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))

##png(paste0("Weight_effect_plots", today, ".png"), height = 7, width = 7, units = "in", res = 400)
jpeg(paste0("Fig4_weight_effect_plots", today, ".jpg"), height = 8.5, width = 7, units = "in", res = 600)
#jpeg(paste0("FigS6_weight_effect_plots_varsnotstd", today, ".jpg"), height = 8.5, width = 7, units = "in", res = 600)
#grid.arrange(p_Year, p_DOY, p_watertemp, p_pressure, p_wind, p_clouds, p_site, ncol = 1)
grid.arrange(p_watertemp, p_pressure, p_wind, p_waterlev, p_clouds, p_moonlight, ncol = 1) #p_watertemp, p_moonlight having issues with the jittering
dev.off()


##----------------------
## AVERAGE YEARLY TRENDS PLOTS
##----------------------
load("weight_year_pred_noenv.RData")  #this is with year as factor and no environmental data
load("weight_year_pred_yrftempd.RData")  #this is with year as factor and year conditioned on watertemp (watertemp as density)
load("weight_year_pred_yrftempc.RData")  #this is with year as factor and year not conditioned on watertemp (watertemp as catchability)
load("wdat.RData")
load("weight_fits_yrftempd.Rdata")

weightenvdat <- subset(weightenvdat, Month != "Oct" & !is.na(wt))

weight_year_pred_noenv$Lake <- factor(weight_year_pred_noenv$Lake, levels = lakes)

weightenvdat$Lake <- factor(weightenvdat$Lake, levels = lakes)

weightenvdat$Lake <- fct_relevel(weightenvdat$Lake,
                                 "BOH", "Furnace", "Feeagh", "Bunaveela")
weight_year_pred_noenv$Lake    <- fct_relevel(weight_year_pred_noenv$Lake,    levels(weightenvdat$Lake))
weight_year_pred_yrftempd$Lake <- fct_relevel(weight_year_pred_yrftempd$Lake, levels(weightenvdat$Lake))
weight_year_pred_yrftempc$Lake <- fct_relevel(weight_year_pred_yrftempc$Lake, levels(weightenvdat$Lake))

# write.csv(weight_year_pred_noenv, file = "weight_year_pred_noenv.csv", row.names = FALSE)
# write.csv(weight_year_pred_yrftempd, file = "weight_year_pred_yrftempd.csv", row.names = FALSE)
# write.csv(weight_year_pred_yrftempc, file = "weight_year_pred_yrftempc.csv", row.names = FALSE)

all_years <- 1987:2025

weightenvdat$Year_num              <- as.numeric(as.character(weightenvdat$Year))
weight_year_pred_noenv$Year_num    <- as.numeric(as.character(weight_year_pred_noenv$Year))
weight_year_pred_yrftempd$Year_num <- as.numeric(as.character(weight_year_pred_yrftempd$Year))
weight_year_pred_yrftempc$Year_num <- as.numeric(as.character(weight_year_pred_yrftempc$Year))

off_box  <- -0.30   # boxplot + red mean CI
off_grn  <- 0       # green diamond
off_blue <- +0.30   # blue triangles

py0 <-
  ggplot()+
  facet_wrap(~Lake, scales = "free", ncol = 2) +
  # first the raw (nominal) data
  #geom_boxplot(data = weightenvdat, 
  #             aes(x = Year_num + off_box, y = wt/Effort, group = Year)) +
  stat_summary(data = weightenvdat,
               aes(x = Year_num + off_box, y = wt/Effort),
               fun.data = "mean_cl_boot",
               colour = "red",
               size = 0.3) +
  # then the values without environmental variables (green star)
  geom_errorbar(data = weight_year_pred_noenv,
                aes(x = Year_num + off_grn, ymin = lwr, ymax = upr, y = yhat),
                width = 0.1, colour = "green") +
  geom_point(data = weight_year_pred_noenv,
             aes(x = Year_num + off_grn, y = yhat),
             size = 2, shape = 18, colour = "green") +  # green diamond
  # then the values with env variables, conditioned on temp (blue triangles)
  geom_errorbar(data = weight_year_pred_yrftempd,
                aes(x = Year_num + off_blue, ymin = lwr, ymax = upr, y = yhat),
                width = 0.1, colour = "blue") +
  geom_point(data = weight_year_pred_yrftempd,
             aes(x = Year_num + off_blue, y = yhat),
             size = 1.5, shape = 17, colour = "blue") +  # blue triangle
  # then the values with env variables, NOT conditioned on temp (red squares)
  # geom_errorbar(data = weight_year_pred_yrftempc,
  #               aes(x = Year_num + off_pur, ymin = lwr, ymax = upr, y = yhat),
  #               width = 0.12, colour = "purple") +
  # geom_point(data = weight_year_pred_yrftempc,
  #            aes(x = Year_num + off_pur, y = yhat),
  #            size = 2, shape = 15, colour = "purple") +   #purple square
  ## Labels + theme
  # scale_x_continuous(
  #   breaks = seq(1987, 2025, 5),
  #   labels = seq(1987, 2025, 5),
  #   limits = c(min(all_years)-0.5, max(all_years)+0.5)) +
  xlab("Year") +
  ylab("Mass per fyke net (kg)") +
  theme(strip.background = element_rect(fill = "grey"))

jpeg(paste0("Fig7_weight_per_year_plots", today, ".jpg"), height = 4.5, width = 10, units = "in", res = 600)
py0
dev.off()



#PLOT comparing raw data to standardized value (conditioned on temp)
y_scales <- list(
  BOH        = scale_y_continuous(limits = c(0, 20)),
  Furnace   = scale_y_continuous(limits = c(0, 3.5)),
  Feeagh    = scale_y_continuous(limits = c(0, 1.75)),
  Bunaveela = scale_y_continuous(limits = c(0, 0.55))
)

#get the right order of the lakes
lake_levels <- c("BOH", "Furnace", "Feeagh", "Bunaveela")
weightenvdat2 <- weightenvdat %>%
  mutate(Lake = factor(Lake, levels = lake_levels))
year_pred2 <- weight_year_pred_yrftempd %>%
  mutate(Lake = factor(Lake, levels = lake_levels))

py1 <-
  ggplot() +
  facet_wrap(~Lake, ncol = 4, scales = "free_y") +
  geom_boxplot(
    data = weightenvdat2,
    aes(x = Year_num, y = wt/Effort, group = Year)
  ) +
  stat_summary(
    data = weightenvdat2,
    aes(x = Year_num, y = wt/Effort),
    fun.data = "mean_cl_boot",
    colour = "red",
    size = 0.3
  ) +
  geom_ribbon(
    data = year_pred2,
    aes(x = Year_num, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.40,
    fill = "lightblue"
  ) +
  geom_line(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    colour = "blue",
    linewidth = 0.5
  ) +
  geom_point(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    size = 2, shape = 15, colour = "blue"
  ) +
  facetted_pos_scales(
    y = y_scales
  ) +
  xlab("Year") +
  ylab("Mass per fyke net (kg)") +
  theme(strip.background = element_rect(fill = "grey"))


py2 <-
  ggplot() +
  facet_wrap(~Lake, ncol = 4, scales = "free_y") +
  geom_ribbon(
    data = year_pred2,
    aes(x = Year_num, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.40,
    fill = "lightblue"
  ) +
  geom_line(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    colour = "blue",
    linewidth = 0.5
  ) +
  geom_point(
    data = year_pred2,
    aes(x = Year_num, y = yhat),
    size = 2, shape = 15, colour = "blue"
  ) +
  xlab("Year") +
  ylab("Mass per fyke net (kg)") +
  theme(strip.background = element_rect(fill = "grey"))



## percentage change between the start and the end of the time series
weight_percent_df <- NULL

for(lake in lakes){
  f0 <- weight_fits_yrftempd[[lake]]
  sub_dat <- subset(weightenvdat, Lake == lake & Month != "Oct" & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
  #only use data with all environmental variables present so same dataset used in all models
  vars <- c("watertemp", "pressure", "wind", "clouds", "moonlight", "waterlev")
  idx <- apply(sub_dat[,vars], 1, function(x){all(!is.na(x))})
  sub_dat <- sub_dat[idx,]
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
  year_temp <- sub_dat |> 
    group_by(Year) |> 
    summarise(watertemp = mean(watertemp, na.rm = TRUE))
  ##
  pred_df2 <- data.frame(Year = year_temp$Year,
                         DOY = means$DOY,
                         pressure = means$pressure,
                         wind = means$wind,
                         clouds = means$clouds,
                         moonlight = means$moonlight,
                         waterlev = means$waterlev,
                         watertemp = year_temp$watertemp,  # <-- CONDITIONED on actual temps
                         fSite = unique(sub_dat$fSite)[1],
                         Effort = 1)
  ##
  Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fSite)"))
  ## draw from posterior - see: https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
  br <- mvrnorm(10000,coef(f0), vcov(f0))
  lnN <- t(Xp %*% t(br))
  ## percentage decline
  #percent_decline <- 100 * (1 - exp(lnN[, 2] - lnN[, 1]))
  if(lake == "BOH"){
    percent_decline <- 100 * (1 - exp(lnN[, 15] - lnN[, 1]))  #because data for 2025 not available
  }
  if(lake == "Bunaveela"){
    percent_decline <- 100 * (1 - exp(lnN[, 21] - lnN[, 1])) #because data for 2025 not available
  }
  if(lake == "Furnace"){
    percent_decline <- 100 * (1 - exp(lnN[, 19] - lnN[, 1])) #because data for 2025 not available
  }
  if(lake == "Feeagh"){
    percent_decline <- 100 * (1 - exp(lnN[, 22] - lnN[, 1])) #because data for 2025 not available
  }
  df <- data.frame(Lake = lake, pd = percent_decline)
  weight_percent_df <- rbind(weight_percent_df, df)
  rm(df)
}

weight_percent_df$Lake <- factor(weight_percent_df$Lake, levels = lakes)

agg0 <- aggregate(pd ~ Lake, FUN = mean, data = weight_percent_df)
agg1 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.025, data = weight_percent_df)
names(agg1)[2] <- c("lwr")
agg2 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.975, data = weight_percent_df)
names(agg2)[2] <- c("upr")

weight_summary_df <- merge(merge(agg0, agg1), agg2)

write.csv(weight_summary_df, file = "weight_percent_decline.csv", row.names = FALSE)

py3 <- ggplot(weight_percent_df, aes(x = pd)) +
  geom_density(fill = "grey") +
  facet_wrap(~ Lake, nrow = 1, scales = "free_y") +
  ylab("Density") +
  geom_vline(data = weight_summary_df, aes(xintercept = pd)) +
  geom_vline(data = weight_summary_df, aes(xintercept = lwr), lty = 3) +
  geom_vline(data = weight_summary_df, aes(xintercept = upr), lty = 3) +
  xlab("Percentage decline in standardised mass trend 1987-2025") +
  coord_cartesian(xlim = c(-50, 100))

##png(paste0("Weight_percent_decline_plots", today, ".png"), height = 8, width = 10, units = "in", res = 400)
jpeg(paste0("Fig8_weight_percent_decline_plots", today, ".jpg"), height = 7, width = 10, units = "in", res = 600)
grid.arrange(py1, py2, py3, nrow =3)
dev.off()











##### Code for plotting interaction terms (moonlight x clouds)  ##############

###################################################################
##  INTERACTION PLOT SLICES ##############
##################################################################

########  COUNT  #################

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

#create slices
slice_props <- c(0.1, 0.5, 0.85)
#slice_props <- c(0.1, 0.5, 0.9)
slice_vals  <- quantile(countenvdat$moonlight, probs = slice_props)

#create predication grid
# grid for clouds × moonlight_slice × Lake
grid_inter <- expand.grid(
  clouds = seq(min(countenvdat$clouds), max(countenvdat$clouds), length = 200),
  moonlight = slice_vals,
  Lake = lakes
)

grid_inter$moonlight_prop <- factor(slice_props, levels = slice_props)[
  match(grid_inter$moonlight, slice_vals)
]

#predict from the gam
vars <- c("Year","DOY","fSite","trap_depth","trap_gradient","trap_number","fchain","watertemp","pressure","wind","waterlev","clouds","moonlight")
vars <- vars[!vars %in% c("fSite", "fchain", "clouds", "moonlight")]
for (v in vars) {
  # Assign lake-specific means
  lake_means <- tapply(countenvdat[[v]], countenvdat$Lake, mean, na.rm = TRUE)
  grid_inter[[v]] <- lake_means[grid_inter$Lake]
}

#add column Effort where value is always 1 and Survey where always is "Russell"
grid_inter$Effort <- 1
grid_inter$survey <- "Russell"

#now for fSite and fchain
grid_inter$fSite <- NA
grid_inter$fchain <- NA
for (lake in lakes) {
  countenvdat1 <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  
  fSite1 <- countenvdat1 %>%
    dplyr::pull(fSite) %>% 
    unique() %>% 
    .[1]
  
  fchain1 <- countenvdat1 %>%
    dplyr::pull(fchain) %>% 
    unique() %>% 
    .[1]
  
  grid_inter$fSite[grid_inter$Lake==lake] <- as.character(fSite1)
  grid_inter$fchain[grid_inter$Lake==lake] <- as.character(fchain1)
}

#write.csv(grid_inter, file = "grid_inter.csv", row.names = TRUE)

#Predict the GAM interaction effect for each lake
bs_year <- "cr"
bs_other <- "tp"


for (lake in lakes) {
  print(lake)
  ## remove October sampling, which was out of the sampling season 
  sub_dat <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
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
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, k = 5, m = 1, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, k= 5, m = 1, bs = bs_other) +
                         s(trap_gradient, k= 5, m = 1, bs = bs_other) +
                         s(trap_number, k= 5, m = 1, bs = bs_other) +
                         s(fchain, bs = "re") +
                         s(watertemp, m = 1, bs = bs_other) +
                         s(pressure, m = 1, bs = bs_other) +
                         s(wind, m = 1, bs = bs_other) +
                         s(clouds, m = 1, bs = bs_other) +
                         s(moonlight, m = 1, bs = bs_other) +
                         ti(clouds, moonlight, m = 1, bs = bs_other) +
                         s(waterlev, m = 1, bs = bs_other) +
                         #guard +
                         offset(log(Effort)))
  }
  ## fit the model
  f0 <- gam(form,
            select = TRUE,
            method = "REML",
            family = nb(),
            data = sub_dat)
  
  pred <- predict(f0, newdata = grid_inter[grid_inter$Lake == lake, ], se.fit = TRUE, type = "terms")
  term_name <- "ti(clouds,moonlight)"   # match summary(gam_model)
  idx <- grid_inter$Lake == lake
  grid_inter$fit[idx] <- pred$fit[, term_name]
  
  #Confidence bands
  grid_inter$se[idx]    <- pred$se.fit[, term_name]
  grid_inter$lower[idx] <- grid_inter$fit[idx] - 2 * grid_inter$se[idx]
  grid_inter$upper[idx] <- grid_inter$fit[idx] + 2 * grid_inter$se[idx]
}

#to get p-values
load("count_fits_yrstempc.RData")

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(count_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

pv <- pvalue_df[pvalue_df$var == "ti(clouds,moonlight)", ]

#drop values of clouds and moonlight outside observed together per lake, so they are not in the plot
ranges <- countenvdat %>%
  group_by(Lake) %>%
  summarise(
    cmin = min(clouds),
    cmax = max(clouds),
    mmin = min(moonlight),
    mmax = max(moonlight)
  )

grid_inter <- grid_inter %>%
  left_join(ranges, by = "Lake") %>%
  filter(
    clouds >= cmin, clouds <= cmax,
    moonlight >= mmin, moonlight <= mmax
  )

#plot
library(ggplot2)

pv2 <- pv %>%
  mutate(text_col = ifelse(p_value < 0.05, "black", "grey50"))

lake_order <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter$Lake   <- factor(grid_inter$Lake,   levels = lake_order)
countenvdat$Lake  <- factor(countenvdat$Lake,  levels = lake_order)
pv2$Lake          <- factor(pv2$Lake,          levels = lake_order)

ggplot(grid_inter,
       aes(x = clouds, y = fit,
           colour = moonlight_prop,
           group = moonlight_prop)) +
  geom_line(linewidth = 1) +
  # rug showing where clouds data exist
  geom_rug(
    data = countenvdat,
    aes(x = clouds),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.4
  ) +
  facet_wrap(~Lake, nrow = 1) +
  scale_colour_viridis_d(name = "Moonlight slice\n(quantile)") +
  labs(
    x = "Clouds",
    y = "Predicted effect"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90")
  ) +
  # add p-value
  geom_text(
    data = pv2,
    aes(
      x = Inf,
      y = Inf,
      label = p_value2,
      color = I(text_col)
    ),
    hjust = 1.2,
    vjust = 2,
    size = 3,
    inherit.aes = FALSE,
    show.legend = FALSE
  )





########  WEIGHT  #################

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

#create slices
slice_props <- c(0.15, 0.5, 0.85)
#slice_props <- c(0.1, 0.5, 0.9)
slice_vals  <- quantile(weightenvdat$moonlight, probs = slice_props)

#create predication grid
# grid for clouds × moonlight_slice × Lake
grid_inter <- expand.grid(
  clouds = seq(min(weightenvdat$clouds), max(weightenvdat$clouds), length = 200),
  moonlight = slice_vals,
  Lake = lakes
)

grid_inter$moonlight_prop <- factor(slice_props, levels = slice_props)[
  match(grid_inter$moonlight, slice_vals)
]

#predict from the gam
vars <- c("Year","DOY","fSite","watertemp","pressure","wind","waterlev","clouds","moonlight")
vars <- vars[!vars %in% c("fSite", "clouds", "moonlight")]
for (v in vars) {
  # Assign lake-specific means
  lake_means <- tapply(weightenvdat[[v]], weightenvdat$Lake, mean, na.rm = TRUE)
  grid_inter[[v]] <- lake_means[grid_inter$Lake]
}

#add column Effort where value is always 1
grid_inter$Effort <- 1

#now for fSite
grid_inter$fSite <- NA
for (lake in lakes) {
  weightenvdat1 <- subset(weightenvdat, Lake == lake & Month != "Oct" & !is.na(wt))
  
  fSite1 <- weightenvdat1 %>%
    dplyr::pull(fSite) %>% 
    unique() %>% 
    .[1]
  
  grid_inter$fSite[grid_inter$Lake==lake] <- as.character(fSite1)
}

#write.csv(grid_inter, file = "grid_inter.csv", row.names = TRUE)

#Predict the GAM interaction effect for each lake
bs_year <- "cr"
bs_other <- "tp"


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
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, m = 1, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 1, bs = bs_other) +
                         s(pressure, m = 1, bs = bs_other) +
                         s(wind, m = 1, bs = bs_other) +
                         s(clouds, m = 1, bs = bs_other) +
                         s(moonlight, m = 1, bs = bs_other) +
                         ti(clouds, moonlight, m = 1, bs = bs_other) +
                         s(waterlev, m = 1, bs = bs_other) +
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
  
  pred <- predict(f0, newdata = grid_inter[grid_inter$Lake == lake, ], se.fit = TRUE, type = "terms")
  term_name <- "ti(clouds,moonlight)"   # match summary(gam_model)
  idx <- grid_inter$Lake == lake
  grid_inter$fit[idx] <- pred$fit[, term_name]
  
  #Confidence bands
  grid_inter$se[idx]    <- pred$se.fit[, term_name]
  grid_inter$lower[idx] <- grid_inter$fit[idx] - 2 * grid_inter$se[idx]
  grid_inter$upper[idx] <- grid_inter$fit[idx] + 2 * grid_inter$se[idx]
}

#to get p-values
load("weight_fits_yrstempc.RData")

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(weight_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

pv <- pvalue_df[pvalue_df$var == "ti(clouds,moonlight)", ]

#drop values of clouds and moonlight outside observed together per lake, so they are not in the plot
ranges <- weightenvdat %>%
  group_by(Lake) %>%
  summarise(
    cmin = min(clouds),
    cmax = max(clouds),
    mmin = min(moonlight),
    mmax = max(moonlight)
  )

grid_inter <- grid_inter %>%
  left_join(ranges, by = "Lake") %>%
  filter(
    clouds >= cmin, clouds <= cmax,
    moonlight >= mmin, moonlight <= mmax
  )

#plot
library(ggplot2)

pv2 <- pv %>%
  mutate(text_col = ifelse(p_value < 0.05, "black", "grey50"))

lake_order <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter$Lake   <- factor(grid_inter$Lake,   levels = lake_order)
countenvdat$Lake  <- factor(countenvdat$Lake,  levels = lake_order)
pv2$Lake          <- factor(pv2$Lake,          levels = lake_order)

ggplot(grid_inter,
       aes(x = clouds, y = fit,
           colour = moonlight_prop,
           group = moonlight_prop)) +
  geom_line(linewidth = 1) +
  # rug showing where clouds data exist
  geom_rug(
    data = countenvdat,
    aes(x = clouds),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.4
  ) +
  facet_wrap(~Lake, nrow = 1) +
  scale_colour_viridis_d(name = "Moonlight slice\n(quantile)") +
  labs(
    x = "Clouds",
    y = "Predicted effect"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90")
  ) +
  # add p-value
  geom_text(
    data = pv2,
    aes(
      x = Inf,
      y = Inf,
      label = p_value2,
      color = I(text_col)
    ),
    hjust = 1.2,
    vjust = 2,
    size = 3,
    inherit.aes = FALSE,
    show.legend = FALSE
  )









###################################################################
##  INTERACTION PLOT 3-D (not using this in the paper) ##############
##################################################################

######## Eel COUNT  ######################################

#The only values varying in the grid are clouds and moonlight

cloud_seq     <- seq(min(countenvdat$clouds),     max(countenvdat$clouds),     length = 100)
moonlight_seq <- seq(min(countenvdat$moonlight),  max(countenvdat$moonlight),  length = 100)
lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter <- expand.grid(
  clouds = cloud_seq,
  moonlight = moonlight_seq,
  Lake = lakes
)

vars <- c("Year","DOY","fSite","trap_depth","trap_gradient","trap_number","fchain","watertemp","pressure","wind","waterlev","clouds","moonlight")
vars <- vars[!vars %in% c("fSite", "fchain", "clouds", "moonlight")]
for (v in vars) {
    # Assign lake-specific means
    lake_means <- tapply(countenvdat[[v]], countenvdat$Lake, mean, na.rm = TRUE)
    grid_inter[[v]] <- lake_means[grid_inter$Lake]
}

#add column Effort where value is always 1 and Survey where always is "Russell"
grid_inter$Effort <- 1
grid_inter$survey <- "Russell"

#now for fSite and fchain
grid_inter$fSite <- NA
grid_inter$fchain <- NA
for (lake in lakes) {
  countenvdat1 <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))

fSite1 <- countenvdat1 %>%
  dplyr::pull(fSite) %>% 
  unique() %>% 
  .[1]

fchain1 <- countenvdat1 %>%
  dplyr::pull(fchain) %>% 
  unique() %>% 
  .[1]

grid_inter$fSite[grid_inter$Lake==lake] <- as.character(fSite1)
grid_inter$fchain[grid_inter$Lake==lake] <- as.character(fchain1)
}

#write.csv(grid_inter, file = "grid_inter.csv", row.names = TRUE)

#Predict the GAM interaction effect for each lake
bs_year <- "cr"
bs_other <- "tp"


for (lake in lakes) {
  print(lake)
  ## remove October sampling, which was out of the sampling season 
  sub_dat <- subset(countenvdat, Lake == lake & Month != "Oct" & !is.na(count))
  sub_dat <- droplevels(sub_dat)
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
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, k = 5, m = 1, bs = bs_other) +
                         s(fSite, bs = "re", k= 3) +
                         s(trap_depth, k= 5, m = 1, bs = bs_other) +
                         s(trap_gradient, k= 5, m = 1, bs = bs_other) +
                         s(trap_number, k= 5, m = 1, bs = bs_other) +
                         s(fchain, bs = "re") +
                         s(watertemp, m = 1, bs = bs_other) +
                         s(pressure, m = 1, bs = bs_other) +
                         s(wind, m = 1, bs = bs_other) +
                         s(clouds, m = 1, bs = bs_other) +
                         s(moonlight, m = 1, bs = bs_other) +
                         ti(clouds, moonlight, m = 1, bs = bs_other) +
                         s(waterlev, m = 1, bs = bs_other) +
                         #guard +
                         offset(log(Effort)))
  }
  ## fit the model
  f0 <- gam(form,
            select = TRUE,
            method = "REML",
            family = nb(),
            data = sub_dat)

pred <- predict(f0, newdata = grid_inter[grid_inter$Lake == lake, ], se.fit = TRUE, type = "terms")
term_name <- "ti(clouds,moonlight)"   # match summary(gam_model)
idx <- grid_inter$Lake == lake
grid_inter$fit[idx] <- pred$fit[, term_name]

#Confidence bands
grid_inter$se[idx]    <- pred$se.fit[, term_name]
grid_inter$lower[idx] <- grid_inter$fit[idx] - 2 * grid_inter$se[idx]
grid_inter$upper[idx] <- grid_inter$fit[idx] + 2 * grid_inter$se[idx]
}

#to get p-values
load("count_fits_yrstempc.RData")

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(count_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

pv <- pvalue_df[pvalue_df$var == "ti(clouds,moonlight)", ]

#drop values of clouds and moonlight outside observed together per lake, so they are not in the plot
ranges <- countenvdat %>%
  group_by(Lake) %>%
  summarise(
    cmin = min(clouds),
    cmax = max(clouds),
    mmin = min(moonlight),
    mmax = max(moonlight)
  )

grid_inter <- grid_inter %>%
  left_join(ranges, by = "Lake") %>%
  filter(
    clouds >= cmin, clouds <= cmax,
    moonlight >= mmin, moonlight <= mmax
  )

#control order of the lakes in the plots
lake_order <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter$Lake   <- factor(grid_inter$Lake,   levels = lake_order)
countenvdat$Lake  <- factor(countenvdat$Lake,  levels = lake_order)
pv$Lake          <- factor(pv$Lake,          levels = lake_order)



#Faceted 2-D surface plot, but very basic
p_surface <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  geom_tile(aes(fill = fit)) +
  geom_contour(aes(z = fit), colour = "black", alpha = 0.4) +
  facet_wrap(~Lake, nrow = 1) +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma") +
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect",
    #title = "Interaction surface: clouds & moonlight"
  ) +
  theme_minimal(base_size = 12)



#plot with more details
p_interaction <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  #optional contour lines
  geom_contour(aes(z = fit), colour = "white", alpha = 0.6) +
  #rug marks for raw data
  geom_rug(
    data = countenvdat,
    aes(x = clouds, y = moonlight),
    inherit.aes = FALSE,
    sides = "b", alpha = 0.5, position = "jitter"
    #alpha = 0.3
  ) +
  #facet by Lake
  facet_wrap(~Lake, nrow = 1) +
  #fill scale
  scale_fill_viridis_c(option = "C") +
  #labels
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect"
  ) +
  #p-values
  geom_text(
    data = pv,
    aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
    hjust = 1.2, vjust = 2,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "black"),
    guide = "none"   # hides the TRUE/FALSE legend
  ) +
  #scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE")) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0, 1, 0, 1), "lines"),
    legend.position = "right")



#another version of the plot, quite similar
p_interaction <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  
  # Shaded surface (analogous to geom_ribbon shading)
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  
  # Contour lines (analogous to thin black smooth)
  geom_contour(aes(z = fit), colour = "black", size = 0.3, alpha = 0.8) +
  
  #facet by Lake
  facet_wrap(~Lake, nrow = 1) +
  
  # Optional rugs to match your x-axis rugs
  geom_rug(
    data = countenvdat,
    aes(x = clouds, y = moonlight),
    inherit.aes = FALSE,
    sides = "b", alpha = 0.5, position = "jitter"
    #alpha = 0.3
    #sides = "b"   # only bottom; remove if too busy
  ) +
  
  # scale_fill_gradient(
  #   low = "grey100",
  #   high = "grey10",
  #   name = paste("Effect")
  # ) +
  
  scale_fill_viridis_c(option = "magma") +
  
  #p-values
  geom_text(
    data = pv,
    aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
    hjust = 1.2, vjust = 2,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "black"),
    guide = "none"   # hides the TRUE/FALSE legend
  ) +
  
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect"
    #title = "Interaction effect of Clouds × Moonlight"
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.margin = unit(c(0, 1, 0, 1), "lines"),
    legend.position = "right"
  )







######## Eel WEIGHT  ######################################

#The only values varying in the grid are clouds and moonlight

cloud_seq     <- seq(min(weightenvdat$clouds),     max(weightenvdat$clouds),     length = 100)
moonlight_seq <- seq(min(weightenvdat$moonlight),  max(weightenvdat$moonlight),  length = 100)
lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter <- expand.grid(
  clouds = cloud_seq,
  moonlight = moonlight_seq,
  Lake = lakes
)

vars <- c("Year","DOY","fSite","watertemp","pressure","wind","waterlev","clouds","moonlight")
vars <- vars[!vars %in% c("fSite", "clouds", "moonlight")]
for (v in vars) {
  # Assign lake-specific means
  lake_means <- tapply(weightenvdat[[v]], weightenvdat$Lake, mean, na.rm = TRUE)
  grid_inter[[v]] <- lake_means[grid_inter$Lake]
}

#add column Effort where value is always 1
grid_inter$Effort <- 1

#now for fSite
grid_inter$fSite <- NA
for (lake in lakes) {
  weightenvdat1 <- subset(weightenvdat, Lake == lake & Month != "Oct" & !is.na(wt))
  
  fSite1 <- weightenvdat1 %>%
    dplyr::pull(fSite) %>% 
    unique() %>% 
    .[1]
  
  grid_inter$fSite[grid_inter$Lake==lake] <- as.character(fSite1)
}

#write.csv(grid_inter, file = "grid_inter.csv", row.names = TRUE)

#Predict the GAM interaction effect for each lake
bs_year <- "cr"
bs_other <- "tp"


for (lake in lakes) {
  print(lake)
  sub_dat <- subset(weightenvdat, Lake == lake & !is.na(wt))
  sub_dat <- droplevels(sub_dat)
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
                         s(Year, m = 1, bs = bs_year) +
                         s(DOY, m = 1, bs = bs_other, k = 5) +
                         s(fSite, bs = "re") +
                         s(watertemp, m = 1, bs = bs_other) +
                         s(pressure, m = 1, bs = bs_other) +
                         s(wind, m = 1, bs = bs_other) +
                         s(clouds, m = 1, bs = bs_other) +
                         s(moonlight, m = 1, bs = bs_other) +
                         ti(clouds, moonlight, m = 1, bs = bs_other) +
                         s(waterlev, m = 1, bs = bs_other) +
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
  
  pred <- predict(f0, newdata = grid_inter[grid_inter$Lake == lake, ], se.fit = TRUE, type = "terms")
  term_name <- "ti(clouds,moonlight)"   # match summary(gam_model)
  idx <- grid_inter$Lake == lake
  grid_inter$fit[idx] <- pred$fit[, term_name]
  
  #Confidence bands
  grid_inter$se[idx]    <- pred$se.fit[, term_name]
  grid_inter$lower[idx] <- grid_inter$fit[idx] - 2 * grid_inter$se[idx]
  grid_inter$upper[idx] <- grid_inter$fit[idx] + 2 * grid_inter$se[idx]
}

#to get p-values
load("weight_fits_yrstempc.RData")

pvalue_df <- NULL

for(lake in lakes){
  tmp <- summary(weight_fits_yrstempc[[lake]])
  tmp <- tmp$s.table
  colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
  df <- as.data.frame(tmp)
  df$var <- rownames(df)
  rownames(df) <- NULL
  df$Lake <- lake
  pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

pv <- pvalue_df[pvalue_df$var == "ti(clouds,moonlight)", ]

#drop values of clouds and moonlight outside observed together per lake, so they are not in the plot
ranges <- weightenvdat %>%
  group_by(Lake) %>%
  summarise(
    cmin = min(clouds),
    cmax = max(clouds),
    mmin = min(moonlight),
    mmax = max(moonlight)
  )

grid_inter <- grid_inter %>%
  left_join(ranges, by = "Lake") %>%
  filter(
    clouds >= cmin, clouds <= cmax,
    moonlight >= mmin, moonlight <= mmax
  )

#control order of the lakes in the plots
lake_order <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

grid_inter$Lake   <- factor(grid_inter$Lake,   levels = lake_order)
countenvdat$Lake  <- factor(weightenvdat$Lake,  levels = lake_order)
pv$Lake          <- factor(pv$Lake,          levels = lake_order)


#Faceted 2-D surface plot, but very basic
p_surface <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  geom_tile(aes(fill = fit)) +
  geom_contour(aes(z = fit), colour = "black", alpha = 0.4) +
  facet_wrap(~Lake, nrow = 1) +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma") +
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect",
    #title = "Interaction surface: ti(clouds, moonlight)"
  ) +
  theme_minimal(base_size = 12)



#plot with more details
p_interaction <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  #optional contour lines
  geom_contour(aes(z = fit), colour = "white", alpha = 0.6) +
  #rug marks for raw data
  geom_rug(
    data = weightenvdat,
    aes(x = clouds, y = moonlight),
    inherit.aes = FALSE,
    sides = "b", alpha = 0.5, position = "jitter"
    #alpha = 0.3
  ) +
  #facet by Lake
  facet_wrap(~Lake, nrow = 1) +
  #fill scale
  scale_fill_viridis_c(option = "C") +
  #labels
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect"
  ) +
  #p-values
  geom_text(
    data = pv,
    aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
    hjust = 1.2, vjust = 2,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "black"),
    guide = "none"   # hides the TRUE/FALSE legend
  ) +
  #scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE")) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0, 1, 0, 1), "lines"),
    legend.position = "right")



#another version of the plot, quite similar
p_interaction <- ggplot(grid_inter, aes(x = clouds, y = moonlight)) +
  
  # Shaded surface (analogous to geom_ribbon shading)
  geom_raster(aes(fill = fit), interpolate = TRUE) +
  
  # Contour lines (analogous to thin black smooth)
  geom_contour(aes(z = fit), colour = "black", size = 0.3, alpha = 0.8) +
  
  #facet by Lake
  facet_wrap(~Lake, nrow = 1) +
  
  # Optional rugs to match your x-axis rugs
  geom_rug(
    data = weightenvdat,
    aes(x = clouds, y = moonlight),
    inherit.aes = FALSE,
    sides = "b", alpha = 0.5, position = "jitter"
    #alpha = 0.3
    #sides = "b"   # only bottom; remove if too busy
  ) +
  
  # scale_fill_gradient(
  #   low = "grey100",
  #   high = "grey10",
  #   name = paste("Effect")
  # ) +
  
  scale_fill_viridis_c(option = "magma") +
  
  #p-values
  geom_text(
    data = pv,
    aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
    hjust = 1.2, vjust = 2,
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c("FALSE" = "darkgrey", "TRUE" = "black"),
    guide = "none"   # hides the TRUE/FALSE legend
  ) +
  
  labs(
    x = "Clouds",
    y = "Moonlight",
    fill = "Effect"
    #title = "Interaction effect of Clouds × Moonlight"
  ) +
  
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.margin = unit(c(0, 1, 0, 1), "lines"),
    legend.position = "right"
  )
