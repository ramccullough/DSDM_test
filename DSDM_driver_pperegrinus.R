# =============================================================================
# DSDM for P.peregrinus
# =============================================================================
# X_train
# beta_surv <- normal(0, 1, dim = K)
# logit_survival_offset <- normal(0, 1)
# 
# X_surv_locs
# 
# # for the survival data
# logit_adult_survival_surv_locs <- X_surv_locs %*% beta_surv
# adult_survival_surv_locs <- ilogit(logit_adult_survival_surv_locs)
# juvenile_survival_surv_locs <- ilogit(logit_adult_survival_surv_locs - logit_survival_offset)
# 
# distribution(survived_adult) <- binomial(trials_adult, adult_survival_surv_locs)
# distribution(survived_juvenile) <- binomial(trials_juvenile, juvenile_survival_surv_locs)
# 
# # survived_adult is a vector of length two (for the two sites) of the combined
# # number of adults that did survive a year, trials_adult is a vector fo length
# # two of the number that were observed (survived + died)
# 
# # for use in the main model
# logit_adult_survival <- X_train %*% beta_surv
# adult_survival <- ilogit(logit_adult_survival)
# juvenile_survival <- ilogit(logit_adult_survival - logit_survival_offset)
# 
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
clean_pp<-function(file, my_pts){
  pp<-read.csv(file)
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
e <- new("Extent", xmin = 143.9, xmax = 145.5, 
         ymin = -38.8, ymax = -37.7)
covariates_cropped <- crop(covariates, e)

# =============================================================================
# Section 3: Set up model and functions
# =============================================================================
# 3.1 choose covariates for survival and fecundity ----------------------------
survival_covs <- c('meananntemp','roads_prox','fire','pop','twi')
fecundity_covs <- c('meananntemp', 'annprecip', 'slope', 'water_prox', 'ndvi')
all_covs <- c(survival_covs, fecundity_covs)
ncov_fecundity <- length(fecundity_covs)
ncov_survival <- length(survival_covs)

# 3.2 Include data from published papers --------------------------------------
default_logit_survival_adult <- qlogis(0.64)
default_logit_survival_juvenile <- qlogis(0.33)
default_log_fecundity <- log(1.23)

# 3.3 Extract covariates at survey locations ----------------------------------
vals <- raster::extract(covariates_cropped, pp_clean[, c('Longitude', 'Latitude')], cellnumbers = TRUE)
dat <- cbind(pp_clean, vals)

# 3.4 remove rows with NA values and scale ------------------------------------
dat_clean <- dat %>%
  na.omit() %>%
  #mutate_at(vars(all_covs), scale) %>%
  mutate_at(vars(all_covs), as.numeric) %>%
  rename(cell_number = cells)

# 3.5 Define priors for model -------------------------------------------------
beta_fecundity <- normal(0, 1/ncov_fecundity, dim = ncov_fecundity)
beta_survival <- normal(0, 1/ncov_survival, dim = ncov_survival)
likelihood_intercept <- variable()

# 3.6 Function for calculating survival, given covariates ---------------------
get_survival <- function(covs){
  x_survival <- covs[, survival_covs]
  survival_logit_adult <- default_logit_survival_adult + x_survival %*% beta_survival
  survival_logit_juvenile <- default_logit_survival_juvenile + x_survival %*% beta_survival
  adult_survival <- ilogit(survival_logit_adult)
  juvenile_survival <- ilogit(survival_logit_juvenile)
  survival <- list(adult = adult_survival, juvenile = juvenile_survival)
}

# 3.7 Function for calculating fecundity, given covariates --------------------
get_fecundity <- function(covs){
  x_fecundity <- covs[, fecundity_covs]
  fecundity_log = default_log_fecundity + x_fecundity %*% beta_fecundity
  fecundity = exp(fecundity_log)
}

# 3.8 Function for calculating lambda, given covariates -----------------------
get_lambda <- function(survival, fecundity){
  top_row <- cbind(fecundity * survival$juvenile,
                   fecundity * survival$adult)
  bottom_row <- cbind(survival$juvenile,
                      survival$adult)
  matrices <- abind(top_row, bottom_row, along = 3)
  iterated <- greta.dynamics::iterate_matrix(matrices, niter = 20)
  lambda <- iterated$lambda
}

# 3.9 Function for calculating probability of presence ------------------------
get_prob <- function(lambda){
  icloglog(likelihood_intercept + log(lambda))
}

# 3.10 Function for calculating count -----------------------------------------
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
prob <- get_prob(lambdas)
rate <- get_rate(lambdas)

# 4.2 Define likelihood function for presence/absence data --------------------
presence <- dat_clean %>%
  mutate_at('Count', as.numeric) %>%
  mutate(PA = ifelse(Count == 0, 0, 1))

presence <- presence[, 'PA']

distribution(presence) <- bernoulli(prob)

# 4.3 Define likelihood function for count data  ------------------------------
abundance <- dat_clean[, 'Count']
distribution(abundance) <- poisson(rate)

# 4.4 Fit model ---------------------------------------------------------------
m <- model(beta_fecundity,
           beta_survival,
           likelihood_intercept)

# 4.5 Run MCMC chain ----------------------------------------------------------
draws <- mcmc(m, n_samples = 500, chains = 2)
# summary(draws)
# mcmc_hist(draws)

# =============================================================================
# Section 5: Get prediction arrays
# =============================================================================
# 5.1 Function to scale covariates
scale_covs <- function (covs, means, sds) {
  cols_scale <- match(names(means), colnames(covs))
  covs_sub <- covs[, cols_scale]
  covs_sub <- sweep(covs_sub, 2, means, "-")
  covs_sub <- sweep(covs_sub, 2, sds, "/")
  covs[, cols_scale] <- covs_sub
  covs
}

# 5.1 Scale bay region covariates by survey site covariates -------------------
cols <- select(dat_clean, all_covs)
dat_means <- colMeans(cols)
dat_sds <- apply(cols, 2, sd)
covariates_predict <- scale_covs(raster::extract(covariates_cropped, seq(1, ncell(covariates_cropped), 1), df = TRUE), dat_means, dat_sds) %>%
  na.omit()

# 5.2 Get prediction arrays ---------------------------------------------------
survival_predict <- get_survival(covariates_predict)
fecundity_predict <- get_fecundity(covariates_predict)
lambdas_predict <- get_lambda(survival_predict, fecundity_predict)
rate_predict <- get_rate(lambdas_predict)
#prob_presence_predict <- 

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
#prob_presence_map <- 

# =============================================================================
# Section 6: Plot prediction maps
# =============================================================================



