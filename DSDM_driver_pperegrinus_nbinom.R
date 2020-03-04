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
library(RColorBrewer)

# 1.2 Set working directory for downloading data  -----------------------------
setwd("C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational/")
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
#r_res <- .0025
#r_ext <- c(141, 150, -39, -33)
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

covariates <- scale(stack(slope,
                          water_prox,
                          fire,
                          # veg_density,
                          ndvi,
                          temp,
                          prec,
                          twi, 
                          pop,
                          roads_prox))

# 2.5 Mask any missing values in covariate stack ------------------------------
missing <- calc(covariates, anyNA)
mask <- covariates[[1]]*0
mask[missing] <- NA
covariates <- mask(covariates, mask)

# 2.6 Reduce size of cells for computation (change to higher res later) -------
covariates <- aggregate(covariates, 10)

# 2.7 Crop stack to Melbourne bay region --------------------------------------
e <- new("Extent", xmin = 143.9, xmax = 145.8, 
         ymin = -38.8, ymax = -37.7)
covariates_cropped <- crop(covariates, e)

# =============================================================================
# Section 3: Set up model and functions
# =============================================================================
# 3.1 choose covariates for survival and fecundity ----------------------------

survival_covs <- c('meananntemp','roads_prox','fire','pop','twi', 'annprecip', 'slope', 'water_prox', 'ndvi')
fecundity_covs <- c('meananntemp','roads_prox','fire','pop','twi', 'annprecip', 'slope', 'water_prox', 'ndvi')
all_covs <- c(survival_covs, fecundity_covs)
ncov_fecundity <- length(fecundity_covs)
ncov_survival <- length(survival_covs)

# 3.2 Extract covariate data from paper sites (Pahl 1987) ---------------------
# Lysterfield: Lat = 37 58'S Long = 145 37'E
# Sandy point: Lat = 38 24'S Long = 145 14'E
Sites <- c('Lysterfield', 'Sandy_Point')
Longitude <- c(145.6167, 145.212436)
Latitude <- c(-37.9667, -38.401404)

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

# 3.6 remove rows with NA values ----------------------------------------------
dat_clean <- dat %>%
  na.omit() %>%
  mutate_at(vars(all_covs), as.numeric) %>%
  rename(cell_number = cells)

# 3.7 Define priors for model -------------------------------------------------
# 3.7a original priors
#beta_fecundity <- normal(0, 1/ncov_fecundity, dim = ncov_fecundity + 1)
#beta_survival <- normal(0, 1/ncov_survival, dim = ncov_survival + 1)
#likelihood_intercept <- variable()

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


# 3.7c Set means and sds for priors -------------------------------------------
# Assume mean + 2*sd = 5
fecundity_prior <- lognormal_prior(1.23, 1.885)
# Assume mean + 2*sd = 0.9
survival_prior <- logitnormal_prior(0.64, 0.13)

dat_survival <- as.matrix(dat_clean[, c(survival_covs, 'intercept')])
dat_fecundity <- as.matrix(dat_clean[, c(fecundity_covs, 'intercept')])

beta_survival_prior <- beta_prior(dat_survival, survival_prior$mean, survival_prior$sd)
beta_fecundity_prior <- beta_prior(dat_fecundity, fecundity_prior$mean, fecundity_prior$sd)

beta_survival <- beta_parameter(beta_survival_prior)
beta_fecundity <- beta_parameter(beta_fecundity_prior)
likelihood_intercept <- variable()

# 3.8 Include survival data from published papers into model ------------------
logit_survival_adult_sites <- x_data[, c(survival_covs, 'intercept')] %*% beta_survival
logit_survival_offset <- variable(lower = 0, upper = Inf) 
survival_adult_sites <- ilogit(logit_survival_adult_sites)
survival_juvenile_sites <- ilogit(logit_survival_adult_sites - logit_survival_offset)

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
  juvenile_survival <- ilogit(survival_logit_adult - logit_survival_offset)
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

# 3.13 Function for calculating probability of presence ------------------------
# get_prob <- function(lambda){
#   icloglog(likelihood_intercept + log(lambda))
# }

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

# 4.2 Define likelihood function for presence/absence data --------------------
# presence <- dat_clean %>%
#   mutate_at('Count', as.numeric) %>%
#   mutate(PA = ifelse(Count == 0, 0, 1))
# 
# presence <- presence[, 'PA']
# 
# distribution(presence) <- bernoulli(prob)

# 4.3 Define likelihood function for count data  ------------------------------
abundance <- dat_clean[, 'Count']

# set prior for dispersion parameter (move to prior section 3.7 when done)
# size <- normal(0, 10, truncation = c(0, Inf))
gamma <- normal(0, 10, truncation = c(0, Inf))
size <- 1/gamma^2

# define probability for negative binomial
nb_prob  <- size / (size + rate)

# negative binomial observation model
distribution(abundance) <- negative_binomial(size, nb_prob)

# 4.4 Fit model ---------------------------------------------------------------
m <- model(beta_fecundity,
           beta_survival,
           likelihood_intercept,
           size)

# 4.5 Run MCMC chain ----------------------------------------------------------
draws <- mcmc(m, n_samples = 500, chains = 4) # increase n chains

# 4.6 Plot summary stats ------------------------------------------------------
setwd("C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Plots/DSDM Outputs/Negative Binomial")
summary(draws)
coda::effectiveSize(draws)
coda::gelman.diag(draws)
mcmc_trace(draws)
mcmc_hist(draws)
mcmc_acf_bar(draws)
mcmc_intervals(draws)

# 4.7 Correlation matrix ------------------------------------------------------
post_cor <- cor(as.matrix(draws), method = 'spearman')
install.packages('corrplot')
library('corrplot')
library(RColorBrewer)
Blues <- colorRampPalette(brewer.pal(9, 'Blues'))((20))
Reds <- colorRampPalette(brewer.pal(9, 'Reds'))(20)

heatmap(x = post_cor, col = Blues, symm = TRUE)

# 4.8 Find model predictions at survey sites for model verification (nbinom) --
nb_prob_draws <- calculate(nb_prob, draws)
nb_prob_matrix <- as.matrix(nb_prob_draws)

size_draws <- calculate(size, draws)
size_matrix <- as.matrix(size_draws)

abundance_sim <- array(NA, dim = dim(nb_prob_matrix))
abundance_sim[] <- rnbinom(length(nb_prob_matrix), size_matrix[], nb_prob_matrix[])

# histogram of simulated abundance
hist(colMeans(abundance_sim), xlab = 'mean simulated abundance', main = 'Simulated abundance, negative binomial')

# 4.9 Compare simulated vs observed data --------------------------------------
ci <- apply(abundance_sim, 2, quantile, c(0.025, 0.975))
points <- cbind(ci[1, ], colMeans(abundance_sim), ci[2,])
plot(abundance, points[, 2], 
     ylim = c(0, 13),
     xlab = "Observed abundance",
     ylab = "Simulated abundance, 95% CI",
     main = "Negative binomial observation model")
arrows(abundance, points[, 1], abundance, points[, 3], length = 0.05, angle = 90, code = 3)
abline(0, 1)

# 4.10 Calculate quantile residuals -------------------------------------------
# package
dharma_obj <- createDHARMa(t(abundance_sim), abundance, integerResponse = TRUE)
hist(residuals(dharma_obj)) 

# manual jitter
abundance_sim_jitter <- jitter(abundance_sim)
abundance_jitter <- jitter(abundance)
p <- colMeans(sweep(abundance_sim_jitter, 2, STATS = abundance_jitter, FUN = `<`))
hist(p, ylim = c(0, 100), xlab = "Randomised quantile residual", main = "Negative Binomial observation model")

# abundance_sim_jitter <- abundance_sim
# abundance_sim_jitter[] <- abundance_sim[] + runif(length(abundance_sim), -0.1, 0.1)
# abundance_jitter <- abundance + runif(length(abundance), -0.1, 0.1)
# p <- colMeans(sweep(abundance_sim_jitter, 2, STATS = abundance_jitter, FUN = `<`))
# hist(p, ylim = c(0, 100), xlab = "Randomised quantile residual", main = "Negative Binomial observation model")

# 4.11 PPC density ------------------------------------------------------------
ppc_dens_overlay(abundance, abundance_sim)

# discretised PPC plot 

# 4.12 Z-scores ---------------------------------------------------------------
# Plot Z score at survey sites against covariate values (vals) at those locations
covs_plot <- dat_clean[, survival_covs]
z_scores <- qnorm(p)

plot(covs_plot[, 1], z_scores, xlab = colnames(covs_plot[1]))
plot(covs_plot[, 2], z_scores, xlab = colnames(covs_plot[2]))
plot(covs_plot[, 3], z_scores, xlab = colnames(covs_plot[3]))
plot(covs_plot[, 4], z_scores, xlab = colnames(covs_plot[4]))
plot(covs_plot[, 5], z_scores, xlab = colnames(covs_plot[5]))
plot(covs_plot[, 6], z_scores, xlab = colnames(covs_plot[6]))
plot(covs_plot[, 7], z_scores, xlab = colnames(covs_plot[7]))
plot(covs_plot[, 8], z_scores, xlab = colnames(covs_plot[8]))
plot(covs_plot[, 9], z_scores, xlab = colnames(covs_plot[9]))

# =============================================================================
# Section 5: Get prediction arrays
# =============================================================================
# 5.1 Get prediction covariates -----------------------------------------------
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
# prob_presence_map <- map_variable(prob_presence_predict, draws) 
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

# 6.2 Plot survival -----------------------------------------------------------
#png(paste("survival", date, sep = "_"))
plot(surv_map, 
     zlim = c(0, 1), 
     main = "Adult survival rate", 
     col = surv_cols,
     axes = FALSE,
     box = FALSE)
#dev.off()

# 6.3 Plot fecundity ----------------------------------------------------------
#png(paste("fecundity", date, sep = "_"))
plot(fec_map, 
     #zlim = c(0, maxValue(fec_map)), 
     main = "Fecundity", 
     col = fec_cols,
     axes = FALSE,
     box = FALSE)
#dev.off()

# 6.4 Plot rate ---------------------------------------------------------------
#png(paste("rate", date, sep = "_"))
plot(rate_map, 
     col = rate_cols,
     main = "Abundance",
     axes = FALSE,
     box = FALSE)
#dev.off()

# 6.5 Plot range --------------------------------------------------------------
#png(paste("range", date, sep = "_"))
range <- lambda_map > 1
plot(range,
     col = c(grey(0.9), 'lightblue'),
     legend = FALSE,
     main = "Range",
     axes = FALSE,
     box = FALSE)
#dev.off()

