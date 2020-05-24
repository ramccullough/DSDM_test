# =============================================================================
# DSDM for P.peregrinus
# =============================================================================
# Section 1: Set up and preliminaries
# Section 2: Download and clean datasets
# Section 3: Set up model priors and functions
# Section 4: Compile model and run mcmc
# Section 5: Get prediction arrays
# Section 6: Plot prediction maps

# =============================================================================
# Section 1: Set up and preliminaries
# =============================================================================
# 1.1 Install packages --------------------------------------------------------
#Install required packages
#install.packages("zoon")
#install.packages("devtools")
#install.packages("bayesplot")

library(raster)
library(zoon)
library(greta)
library(bayesplot)
library(dplyr)
library(rgeos)
library(tidyr)
library('corrplot')
library(RColorBrewer)

set.seed(20200420)

# 1.2 Set working directory for downloading data  -----------------------------
setwd("~/")
#setwd("C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational/")
rel_path_input <- 'pp_grids/pp_grids/'

# =============================================================================
# Section 2: Download and clean datasets
# =============================================================================
# 2.1 Function for preparing survey site coordinate data ----------------------
prepare_coords_data <- function (coordinates_file) {
  
  attr <- rgdal::readOGR(coordinates_file, 'Survey_points') %>%
    attributes()
  
  coords <- cbind(attr$coords[, 1], attr$coords[, 2], attr$data)[, c(1:3)] 
  colnames(coords) <- c('Longitude', 'Latitude', 'Site_point')
  coords$Site_point <- as.character(coords$Site_point)
  
  my_pts <- tidyr::separate(coords, 'Site_point', c('Site', 'Point'))
  my_pts$Site <- as.integer(my_pts$Site)
  return(my_pts)
  
}

# 2.2 Function for cleaning survey data into required format ------------------
clean_pp <- function(file, my_pts){
  pp <- read.csv(file)
  pp_clean <- pp %>%
    mutate(RT = Species == "RT") %>%
    group_by(Site, Point) %>%
    summarise(Count = sum(RT)) %>% 
    merge(my_pts, by = c('Site', 'Point'))
  return(pp_clean)
}

# 2.3 Load site and survey data -----------------------------------------------
my_pts <- prepare_coords_data(coordinates_file = 'Survey_points.kml')
pp_clean <- clean_pp('survey_data.csv', my_pts)

# 2.3 Load environmental covariates -------------------------------------------
r <- aggregate(shapefile('pp_grids/pp_grids/cma.shp'), dissolve = TRUE) %>%
  as("SpatialPolygonsDataFrame")

# 2.4 Create raster stack of environmental covariates -------------------------
slope <- mask(raster(x = paste0(rel_path_input, 'slope.tif')), r)
water_prox <- mask(raster(x = paste0(rel_path_input, 'water_prox.tif')), r)
fire <- mask(raster(x = paste0(rel_path_input, 'fire.tif')), r)
# veg_density <- mask(raster(x = paste0(rel_path_input, 'sveg.tif')), r)
ndvi <- mask(raster(x = paste0(rel_path_input, 'ndvi.tif')), r) #normalised vegetation 
temp <- mask(raster(x = paste0(rel_path_input, 'meananntemp.tif')), r)
prec <- mask(raster(x = paste0(rel_path_input, 'annprecip.tif')), r)
twi <- mask(raster(x = paste0(rel_path_input, 'twi.tif')), r) #terrain wetness
pop <- mask(raster(x = paste0(rel_path_input, 'pop.tif')), r) #human population
roads_prox <- mask(raster(x = paste0(rel_path_input, 'roads_prox.tif')), r)

#setwd('C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational/gl_grumpv1_urextent_bil_30')
setwd("~/urban_layer")
global_urban <- raster('glurextents.bil')
urban <- resample(global_urban, slope)

#setwd('C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational/vic_temp_data')
setwd("~/temperature_layer")
vic_temp_days <- raster('av_days_data_10y.tif')
temp_days <- resample(vic_temp_days, slope)

roads_sqrt <- sqrt(roads_prox)

covariates <- scale(stack(slope,
                          water_prox,
                          fire,
                          # veg_density,
                          ndvi,
                          temp,
                          prec,
                          twi, 
                          pop,
                          #roads_prox,
                          roads_sqrt,
                          urban,
                          temp_days))

# 2.5 Mask any missing values in covariate stack ------------------------------
missing <- calc(covariates, anyNA)
mask <- covariates[[1]]*0
mask[missing] <- NA
covariates <- mask(covariates, mask)

# 2.6 Reduce size of cells for computation (change to higher res later) -------
# covariates <- aggregate(covariates, 10)

# 2.7 Crop stack to Melbourne bay region --------------------------------------
e <- new("Extent", xmin = 143.9, xmax = 145.8, 
         ymin = -38.8, ymax = -37.7)
covariates_cropped <- crop(covariates, e)

names(covariates_cropped[[9]]) <- 'roads_sqrt'
names(covariates_cropped[[10]]) <- 'urban'
names(covariates_cropped[[11]]) <- 'temp_days'

# =============================================================================
# Section 3: Set up model and functions
# =============================================================================
# 3.1 choose covariates for survival and fecundity ----------------------------
survival_covs <- c('meananntemp','roads_sqrt','fire','pop','twi', 'annprecip', 'slope', 'water_prox', 'ndvi', 'urban', 'temp_days')
fecundity_covs <- c('meananntemp','roads_sqrt','fire','pop','twi', 'annprecip', 'slope', 'water_prox', 'ndvi', 'urban', 'temp_days')
all_covs <- c(survival_covs, fecundity_covs)
ncov_fecundity <- length(fecundity_covs)
ncov_survival <- length(survival_covs)

# 3.2 Extract covariate data from paper sites (Pahl 1987) ---------------------
# Lysterfield: Lat = 37 58'S Long = 145 37'E
# Sandy point: Lat = 38 24'S Long = 145 14'E
Sites <- c('Lysterfield', 'Sandy_Point')
Longitude <- c(145.6167, 145.212436)
Latitude <- c(-37.9667, -38.401404)

# -38.408055, 145.23998

site_locs <- data.frame(Sites, Longitude, Latitude)
x_data <- raster::extract(covariates, site_locs[, c('Longitude', 'Latitude')], cellnumbers = TRUE)
x_data <- as.data.frame(x_data) %>%
  cbind(intercept = 1)

# 3.3 Include survival data from paper sites (Pahl 1987) ----------------------
adult_trials <- c(54, 65)
adult_survived <- c(23, 34)
juvenile_trials <- c(98, 82)
juvenile_survived <- c(31, 31)

survival_adult <- data.frame(Sites, adult_trials, adult_survived)
survival_juvenile <- data.frame(Sites, juvenile_trials, juvenile_survived)

# 3.4 Include fecundity data from paper sites (Pahl 1987) ---------------------
breeding_females <- c(74, 107)
offspring <- c(176, 239)
female_offspring <- c(88, 120)
fecundity_sites <- data.frame(Sites, breeding_females, offspring, female_offspring)

# 3.5 Extract covariates at survey locations ----------------------------------
vals <- raster::extract(covariates_cropped, pp_clean[, c('Longitude', 'Latitude')], cellnumbers = TRUE)
dat <- cbind(pp_clean, vals, intercept = 1)

# 3.6 remove rows with NA values ?(and scale) ---------------------------------
dat_clean <- dat %>%
  na.omit() %>%
  mutate_at(vars(all_covs), as.numeric) %>%
  rename(cell_number = cells)

# 3.7 Define priors for model -------------------------------------------------
# 3.7a original priors
# beta_fecundity <- normal(0, 1/ncov_fecundity, dim = ncov_fecundity + 1)
# beta_survival <- normal(0, 1/ncov_survival, dim = ncov_survival + 1)
# likelihood_intercept <- variable()

# 3.7b MVN prior
beta_prior <- function(X, mean, sd) {
  
  # use weights to increase variance for each datapoint to spread
  # prior over multiple datapoints, without adding too much extra data
  
  # convert the means and sds into MVN parameters
  n <- nrow(X)
  
  if (length(mean) == 1) {
    mean <- rep(mean, n)
  }
  if (length(sd) == 1) {
    sd <- rep(sd, n)
  }
  
  var <- sd ^ 2
  cov <- diag(var)
  
  # invert X
  iX <- MASS::ginv(X)
  
  # solve for beta
  list(mean = iX %*% mean,
       Sigma = iX %*% cov %*% t(iX))
  
}

beta_parameter <- function(beta_prior) {
  L <- t(chol(beta_prior$Sigma))
  beta_raw <- normal(0, 1, dim = length(beta_prior$mean))
  beta_prior$mean + L %*% beta_raw
}

lognormal_prior <- function(mean, sd) {
  
  var <- sd ^ 2
  list(
    mean = log((mean ^ 2) / sqrt(var + mean ^ 2)),
    sd = sqrt(log(1 + var / (mean ^ 2)))
  )
  
}

logitnormal_prior <- function(mean, sd) {
  
  # given the parameters of the distribution, estimate the moments and hence the
  # mean and variance via the quasi Monte Carlo method listed on wikipedia
  logitnormal_mean_var <- function (mu, sd, K = 1000) {
    i <- 1:(K - 1)
    p <- plogis(qnorm(i/K, mu, sd))
    m1 <- mean(p)
    m2 <- mean(p ^ 2)
    c(mean = m1, variance = m2 - m1 ^ 2)
  }
  
  # compute the root mean squared error between the logitnormal mean and
  # variance from this (transformed) set of parameters, and the observed mean
  # and variance
  obj <- function (par, observed) {
    mu <- par[1]
    sd <- exp(par[2])
    expected <- logitnormal_mean_var(mu, sd)
    mean((expected - observed) ^ 2)
  }
  
  # optimise this to find the parameter set that minimised the mean squared error
  op <- optim(par = c(0, 0), obj, observed = c(mean, sd ^ 2))
  
  # return the mean and variance on the logit scale
  list(mean = op$par[1],
       sd = exp(op$par[2]))
  
}

delta_prior <- function (logit_S_A,
                         logit_S_J,
                         n = 1e5) {
  
  list(mean = logit_S_A$mean - logit_S_J$mean,
       sd = sqrt(logit_S_A$sd ^ 2 + logit_S_J$sd ^ 2))
}


# 3.7c Set means and sds for priors -------------------------------------------
fecundity_prior <- lognormal_prior(1.23, 1)
survival_prior <- logitnormal_prior(0.64, 0.1)
survival_juvenile_prior <- logitnormal_prior(0.33, 0.1)

dat_survival <- as.matrix(dat_clean[, c(survival_covs, 'intercept')])
dat_fecundity <- as.matrix(dat_clean[, c(fecundity_covs, 'intercept')])

beta_survival_prior <- beta_prior(dat_survival, survival_prior$mean, survival_prior$sd)
beta_fecundity_prior <- beta_prior(dat_fecundity, fecundity_prior$mean, fecundity_prior$sd)
delta_prior <- delta_prior(survival_prior, survival_juvenile_prior)

beta_survival <- beta_parameter(beta_survival_prior)
beta_fecundity <- beta_parameter(beta_fecundity_prior)
delta <- normal(delta_prior$mean, delta_prior$sd)
likelihood_intercept <- variable()

# 3.8 Include survival data from published papers into model ------------------
logit_survival_adult_sites <- x_data[, c(survival_covs, 'intercept')] %*% beta_survival
survival_adult_sites <- ilogit(logit_survival_adult_sites)
survival_juvenile_sites <- ilogit(logit_survival_adult_sites - delta)

distribution(survival_adult$adult_survived) <- binomial(survival_adult$adult_trials, survival_adult_sites)
distribution(survival_juvenile$juvenile_survived) <- binomial(survival_juvenile$juvenile_trials, survival_juvenile_sites)

# 3.9 Include fecundity data from published papers into model -----------------
log_fecundity = x_data[, c(fecundity_covs, 'intercept')] %*% beta_fecundity
fec_rate <- exp(log_fecundity)

distribution(fecundity_sites$female_offspring) <- poisson(fec_rate * fecundity_sites$breeding_females)

# 3.10 Function for calculating survival, given covariates ---------------------
get_survival <- function(covs){
  x_survival <- covs[, c(survival_covs, 'intercept')]
  survival_logit_adult <- x_survival %*% beta_survival
  adult_survival <- ilogit(survival_logit_adult)
  #juvenile_survival <- ilogit(survival_logit_adult - logit_survival_offset)
  juvenile_survival <- ilogit(survival_logit_adult - delta)
  survival <- list(adult = adult_survival, juvenile = juvenile_survival)
}

# 3.11 Function for calculating fecundity, given covariates --------------------
get_fecundity <- function(covs){
  x_fecundity <- covs[, c(fecundity_covs, 'intercept')]
  fecundity_log = x_fecundity %*% beta_fecundity
  fecundity = exp(fecundity_log)
}

# 3.12 Function for calculating lambda, given covariates -----------------------
get_lambda <- function(survival, fecundity){
  top_row <- cbind(fecundity * survival$juvenile,
                   fecundity * survival$adult)
  bottom_row <- cbind(survival$juvenile,
                      survival$adult)
  matrices <- abind(top_row, bottom_row, along = 3)
  iterated <- greta.dynamics::iterate_matrix(matrices, niter = 20)
  lambda <- iterated$lambda
}

# 3.14 Function for calculating count -----------------------------------------
get_rate <- function(lambda){
  exp(likelihood_intercept)*lambda
}

# =============================================================================
# Section 4: Compile model and run mcmc
# =============================================================================
# 4.1 Get survival and fecundity for survey sites -----------------------------
survival <- get_survival(dat_clean)
fecundity <- get_fecundity(dat_clean)
lambdas <- get_lambda(survival, fecundity)
# prob <- get_prob(lambdas)
rate <- get_rate(lambdas)

# 4.2 Define likelihood function for count data  ------------------------------
abundance <- dat_clean[, 'Count']
distribution(abundance) <- poisson(rate)

# 4.3 Fit model ---------------------------------------------------------------
m <- model(beta_fecundity,
           beta_survival,
           likelihood_intercept,
           delta)

# 4.4 Run MCMC chain ----------------------------------------------------------
draws <- greta::mcmc(m, n_samples = 1000, chains = 4) 

# Set folder for all output plots
# Create new folder for each day
setwd("~/Plots/Poisson/20200521")

outputs <- list(
  r_hat = coda:::gelman.diag(draws, multivariate = FALSE),
  n_eff = coda:::effectiveSize(draws),
  summary = summary(draws))

# 4.5 Plot summary stats ------------------------------------------------------
# Get colour scheme
Blues <- c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")
Reds <- c("#DCBCBC", "#C79999",
          "#B97C7C", "#A25050", "#8F2727", "#7C0000")

# Change variable names for plot titles
varnames(draws) <- c("beta^{F[1]}", 
                     "beta^{F[2]}",
                     "beta^{F[3]}",
                     "beta^{F[4]}",
                     "beta^{F[5]}",
                     "beta^{F[6]}",
                     "beta^{F[7]}",
                     "beta^{F[8]}",
                     "beta^{F[9]}",
                     "beta^{F[10]}",
                     "beta^{F[11]}",
                     "beta^{F[0]}",
                     "beta^{S[1]}", 
                     "beta^{S[2]}", 
                     "beta^{S[3]}",
                     "beta^{S[4]}",
                     "beta^{S[5]}",
                     "beta^{S[6]}",
                     "beta^{S[7]}",
                     "beta^{S[8]}",
                     "beta^{S[9]}",
                     "beta^{S[10]}",
                     "beta^{S[11]}",
                     "beta^{S[0]}",
                     "alpha",
                     "delta")

# Get mcmc data as tbl
mcmc_data <- mcmc_trace_data(draws)

# Trace plot
beta_s_data <- mcmc_data %>%
  filter(grepl('S', parameter))

beta_f_data <-mcmc_data %>%
  filter(grepl('F', parameter))

other_data <- mcmc_data %>%
  filter(!grepl('beta', parameter))

# Trace plots
png("trace_surv_params.png",
    height = 1600,
    width = 1600)
ggplot(data = beta_s_data, aes(x = iteration, y = value)) +
  geom_line(aes(colour = chain)) +
  facet_wrap(~parameter, labeller = label_parsed) +
  scale_colour_manual(name = "", values = c(Blues[[6]], Blues[[4]], Blues[[3]], Blues[[1]])) +
  guides(colour = guide_legend(title = "Chain")) +
  labs(x = "Iteration", y = "Parameter value") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()


png("trace_fec_params.png",
    height = 1600,
    width = 1600)
ggplot(data = beta_f_data, aes(x = iteration, y = value)) +
  geom_line(aes(colour = chain)) +
  facet_wrap(~parameter, labeller = label_parsed) +
  scale_colour_manual(name = "", values = c(Blues[[6]], Blues[[4]], Blues[[3]], Blues[[1]])) +
  guides(colour = guide_legend(title = "Chain")) +
  labs(x = "Iteration", y = "Parameter value") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

png("trace_other_params.png",
    height = 530,
    width = 800)
ggplot(data = other_data, aes(x = iteration, y = value)) +
  geom_line(aes(colour = chain)) +
  facet_wrap(~parameter, labeller = label_parsed) +
  scale_colour_manual(name = "", values = c(Blues[[6]], Blues[[4]], Blues[[3]], Blues[[1]])) +
  guides(colour = guide_legend(title = "Chain")) +
  labs(x = "Iteration", y = "Parameter value") + 
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

# Histograms
png("histogram_surv_params.png",
    height = 1600,
    width = 1600)
ggplot() + 
  geom_histogram(data = beta_s_data, binwidth = 0.01, aes(x = value, y = ..density.., fill = "HMC draws"), colour = Blues[[3]], size = 0.25) +
  scale_fill_manual(name = "", values = Blues[[2]]) +
  facet_wrap(~parameter, labeller = label_parsed, scales = "free_y") +
  guides(fill = guide_legend(order = 1)) +
  labs(x = "Parameter value", y = "Density") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

png("histogram_fec_params.png",
    height = 1600,
    width = 1600)
ggplot() + 
  geom_histogram(data = beta_f_data, binwidth = 0.01, aes(x = value, y = ..density.., fill = "HMC draws"), colour = Blues[[3]], size = 0.25) +
  scale_fill_manual(name = "", values = Blues[[2]]) +
  facet_wrap(~parameter, labeller = label_parsed, scales = "free_y") +
  guides(fill = guide_legend(order = 1)) +
  labs(x = "Parameter value", y = "Density") + 
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

png("histogram_other_params.png",
    height = 530,
    width = 800)
ggplot() + 
  geom_histogram(data = other_data, binwidth = 0.01, aes(x = value, y = ..density.., fill = "HMC draws"), colour = Blues[[3]], size = 0.25) +
  scale_fill_manual(name = "", values = Blues[[2]]) +
  facet_wrap(~parameter, labeller = label_parsed, scales = "free_y") +
  guides(fill = guide_legend(order = 1)) +
  labs(x = "Parameter value", y = "Density") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()


# 4.6 Correlation matrix ------------------------------------------------------
# post_cor <- cor(as.matrix(draws), method = 'spearman')
# install.packages('corrplot')
# library('corrplot')
# library(RColorBrewer)
# Blues <- colorRampPalette(brewer.pal(9, 'Blues'))((20))
# Reds <- colorRampPalette(brewer.pal(9, 'Reds'))(20)
# 
# heatmap(x = post_cor, col = Blues, symm = TRUE)

# 4.7 Find model predictions at survey sites for model verification (pois) ----
rate_draws <- calculate(rate, draws)
rate_matrix <- as.matrix(rate_draws)

abundance_sim <- array(NA, dim = dim(rate_matrix))
abundance_sim[] <- rpois(length(rate_matrix), rate_matrix[])

# (manual)
abundance_sim_jitter <- jitter(abundance_sim)
abundance_jitter <- jitter(abundance)
p <- colMeans(sweep(abundance_sim_jitter, 2, STATS = abundance_jitter, FUN = `>`))

bin_breaks <- seq(0, 1, 0.1)

png("pvalue.png",
    height = 800,
    width = 1200)
p_data <- data.frame(quantile = p)
ggplot(p_data) +
  geom_histogram(breaks = bin_breaks, aes(x = quantile), colour = Blues[[3]], fill = Blues[[2]]) +
  labs(x = "Posterior predictive p-value", y = "Count") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

# QQ plot
z_data <- data.frame(z_score = qnorm(p))

png("qqplot.png",
    height = 800,
    width = 1200)
ggplot(z_data) +
  geom_qq(aes(sample = z_score), colour = Blues[[5]], size = 3) +
  geom_abline() +
  labs(x = "Standard normal quantiles", y = "Model prediction quantile residuals") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()


# 4.8 PPC density ------------------------------------------------------------
ppc_data <- ppc_data(abundance, abundance_sim)

png("ppc_overlay.png",
    height = 800,
    width = 1200)
ggplot(ppc_data) +
  stat_ecdf(data = function(x) dplyr::filter(x, !.data$is_y), aes(x = value,  group = rep_id, colour = "Model predictions"), size = 0.25) +
  stat_ecdf(data = function(x) dplyr::filter(x, .data$is_y), aes(x = value, colour = "Observed data"), size = 1.0) +
  scale_colour_manual(name = "", values = c(Blues[[2]], Blues[[6]])) +
  guides(colour = guide_legend(title = "")) +
  labs(x = "Abundance value", y = "Cumulative density") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

# 4.9 Z-scores ---------------------------------------------------------------
##  Plot Z score at survey sites against covariate values (vals) at those locations
# covs_plot <- dat_clean[, survival_covs]
# z_scores <- qnorm(p)
# 
# # z_scores <- (colMeans(abundance_sim) - mean(colMeans(abundance_sim)))/sd(abundance_sim)
# 
# plot(covs_plot[, 1], z_scores, xlab = colnames(covs_plot[1]))
# plot(covs_plot[, 2], z_scores, xlab = colnames(covs_plot[2]))
# plot(covs_plot[, 3], z_scores, xlab = colnames(covs_plot[3]))
# plot(covs_plot[, 4], z_scores, xlab = colnames(covs_plot[4]))
# plot(covs_plot[, 5], z_scores, xlab = colnames(covs_plot[5]))
# plot(covs_plot[, 6], z_scores, xlab = colnames(covs_plot[6]))
# plot(covs_plot[, 7], z_scores, xlab = colnames(covs_plot[7]))
# plot(covs_plot[, 8], z_scores, xlab = colnames(covs_plot[8]))
# plot(covs_plot[, 9], z_scores, xlab = colnames(covs_plot[9]))
# plot(covs_plot[, 10], z_scores, xlab = colnames(covs_plot[10]))
# plot(covs_plot[, 11], z_scores, xlab = colnames(covs_plot[11]))

# =============================================================================
# Section 5: Get prediction arrays
# =============================================================================
# 5.1 Extract covariates at all prediction sites ------------------------------
covariates_predict <- raster::extract(covariates_cropped, seq(1, ncell(covariates_cropped), 1)) %>%
  cbind(intercept = 1) %>%
  na.omit()

# 5.2 Get prediction arrays ---------------------------------------------------
survival_predict <- get_survival(covariates_predict)
fecundity_predict <- get_fecundity(covariates_predict)
lambdas_predict <- get_lambda(survival_predict, fecundity_predict)
# prob_presence_predict <- get_prob(lambdas_predict)
rate_predict <- get_rate(lambdas_predict)

# 5.3 Get cell IDs of non-NA cells in prediction region -----------------------
idx <- which(!is.na(getValues(covariates_cropped[[1]])))

# 5.4 Function for mapping mean of posterior sample for map grid cells --------
map_variable <- function (greta_array, draws) {
  vals <- calculate(greta_array, draws)
  vals_mat <- as.matrix(vals)
  mean <- colMeans(vals_mat)
  map <- covariates_cropped[[1]]
  map[idx] <- mean
  map
}

# 5.5 Get values for mapping --------------------------------------------------
lambda_map <- map_variable(lambdas_predict, draws)
rate_map <- map_variable(rate_predict, draws)
fec_map <- map_variable(fecundity_predict, draws)
surv_map <- map_variable(survival_predict$adult, draws)

# =============================================================================
# Section 6: Plot prediction maps
# =============================================================================
# 6.1 Set up colours for plotting ---------------------------------------------
YlGnBu <- brewer.pal(9, 'YlGnBu')[c(1, 1, 2, 3, 4, 5, 6, 7, 7)]
YlGn <- brewer.pal(9, 'YlGn')
RdPu <- brewer.pal(9, 'RdPu')
Greens <- brewer.pal(9, 'Greens')

rate_cols <- colorRampPalette(YlGnBu)(1000)
surv_cols <- colorRampPalette(YlGn)(1000)
fec_cols <- colorRampPalette(RdPu)(1000)


# combine all plots into same figure
png('poisson_outputs.png',
    width = 2400,
    height = 1600,
    pointsize = 45)

par(mfrow=c(2,2),
    oma = c(1, 1, 1, 1),
    mai = c(0, 0, 0, 0),
    mar = c(1, 0, 2, 0))

# 6.2 Plot survival -----------------------------------------------------------
plot(surv_map, 
     zlim = c(0, 1), 
     main = "Adult survival rate", 
     col = surv_cols,
     axes = FALSE,
     box = FALSE)

# 6.3 Plot fecundity ----------------------------------------------------------
plot(fec_map, 
     zlim = c(0, 5), 
     main = "Fecundity", 
     col = fec_cols,
     axes = FALSE,
     box = FALSE)

# 6.4 Plot rate ---------------------------------------------------------------
plot(rate_map,
     zlim = c(0, 2.0),
     col = rate_cols,
     main = "Abundance",
     axes = FALSE,
     box = FALSE)

# 6.5 Plot range --------------------------------------------------------------
range <- lambda_map > 1
plot(range,
     col = c(grey(0.9), Blues[6]),
     legend = FALSE,
     main = "Range",
     axes = FALSE,
     box = FALSE)
dev.off()

