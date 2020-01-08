## DSDM for P.peregrinus
#=========================
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

# Set working directory for data
setwd("C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational/")
rel_path_input <- 'pp_grids/pp_grids/'

# Function for coordinates data
prepare_coords_data <- function (coordinates_file) {
  
  attr <- rgdal::readOGR(coordinates_file, 'Survey_points') %>%
    attributes()
  
  coords <- cbind(attr$coords[,1], attr$coords[,2], attr$data)[,c(1:3)] 
  colnames(coords) <- c('Longitude', 'Latitude', 'Site_point')
  coords$Site_point <- as.character(coords$Site_point)
  
  my_pts <- tidyr::separate(coords, 'Site_point', c('Site', 'Point'))
  my_pts$Site <- as.integer(my_pts$Site)
  return(my_pts)
  
}

# Function for cleaning survey data
clean_pp<-function(file, my_pts){
  pp<-read.csv(file)
  pp_clean <- pp %>%
    mutate(RT = Species == "RT") %>%
    group_by(Site, Point) %>%
    summarise(Count = sum(RT)) %>% 
    merge(my_pts, by = c('Site', 'Point'))
  
  return(pp_clean)
}

# Load survey data
my_pts <- prepare_coords_data(coordinates_file = 'Survey_points.kml')
pp_clean <- clean_pp('survey_data.csv', my_pts)

#load environmental covariates
#r_res <- .0025
#r_ext <- c(141, 150, -39, -33)

r <- aggregate(shapefile('pp_grids/pp_grids/cma.shp'), dissolve = TRUE) %>%
  as("SpatialPolygonsDataFrame")

#create raster stack of environmental covariates
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

#stack rasters
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

# Mask any missing values in covariate stack (eg ocean in bay)
missing <- calc(covariates, anyNA)
mask <- covariates[[1]]*0
mask[missing] <- NA
covariates <- mask(covariates, mask)

# reduce size for computation
covariates <- aggregate(covariates, 10)

# Crop to Melbourne bay region
e <- new("Extent", xmin = 143.9, xmax = 145.5, 
         ymin = -38.8, ymax = -37.7)
covariates_cropped <- crop(covariates, e)

##Refine/edit
#choose covariates for survival and fecundity.
survival_covs <- c('meananntemp','roads_prox','fire','pop','twi')
fecundity_covs <- c('meananntemp', 'annprecip', 'slope', 'water_prox', 'ndvi')
all_covs <- c(survival_covs, fecundity_covs)
ncov_fecundity <- length(fecundity_covs)
ncov_survival <- length(survival_covs)

# Enter survival/fecundity data from papers here
# Placeholder numbers for now for testing
default_logit_survival_adult <- qlogis(0.64)
default_logit_survival_juvenile <- qlogis(0.33)
default_log_fecundity <- log(1.23)

# Extract covariates at survey locations
vals <- raster::extract(covariates_cropped, pp_clean[, c('Longitude', 'Latitude')], cellnumbers = TRUE)
dat <- cbind(pp_clean, vals)

#remove rows with NA values for covariates, and scale
dat_clean <- dat %>%
  na.omit() %>%
  #mutate_at(vars(all_covs), scale) %>%
  mutate_at(vars(all_covs), as.numeric) %>%
  rename(cell_number = cells)

beta_fecundity <- normal(0, 1/ncov_fecundity, dim = ncov_fecundity)
beta_survival <- normal(0, 1/ncov_survival, dim = ncov_survival)
likelihood_intercept <- variable()

# calculate survival, given covariates
#survival
get_survival <- function(covs){
  x_survival <- covs[, survival_covs]
  survival_logit_adult <- default_logit_survival_adult + x_survival %*% beta_survival
  survival_logit_juvenile <- default_logit_survival_juvenile + x_survival %*% beta_survival
  adult_survival <- ilogit(survival_logit_adult)
  juvenile_survival <- ilogit(survival_logit_juvenile)
  survival <- list(adult = adult_survival, juvenile = juvenile_survival)
}

#fecundity
get_fecundity <- function(covs){
  x_fecundity <- covs[, fecundity_covs]
  fecundity_log = default_log_fecundity + x_fecundity %*% beta_fecundity
  fecundity = exp(fecundity_log)
}

#lambda
get_lambda <- function(survival, fecundity){
  top_row <- cbind(fecundity * survival$juvenile,
                   fecundity * survival$adult)
  bottom_row <- cbind(survival$juvenile,
                      survival$adult)
  matrices <- abind(top_row, bottom_row, along = 3)
  iterated <- greta.dynamics::iterate_matrix(matrices, niter = 20)
  lambda <- iterated$lambda
}

#probability of presence
get_prob <- function(lambda){
  icloglog(likelihood_intercept + log(lambda))
}

#poisson rate
get_rate <- function(lambda){
  exp(likelihood_intercept)*lambda
}

#get survival and fecundity
survival <- get_survival(dat_clean)
fecundity <- get_fecundity(dat_clean)
lambdas <- get_lambda(survival, fecundity)
prob <- get_prob(lambdas)
rate <- get_rate(lambdas)

#likelihood for PA
# To finish
#presence<- dat_clean %>%
  #mutate_at('Count', as.logical) %>%
  #mutate_at('Count', as.numeric)
#distribution(presence)<-binomial(prob)

#likelihood for counts
abundance <- dat_clean[, 'Count']
distribution(abundance) <- poisson(rate)

#fit model
m <- model(beta_fecundity,
           beta_survival,
           likelihood_intercept)

draws <- mcmc(m, n_samples = 500, chains = 2)
# summary(draws)
# mcmc_hist(draws)

# Predicting from the model

# Get covariate values for all pixels in raster
# Using only cropped Melbourne bay region

# scale prediction covariates by mean and sd of survey site covariates
# Columns of covariates only
# dat_means <- dat %>%
#   na.omit() %>%
#   select(all_covs) %>%
#   summarise_all(mean) %>%
#   as.numeric()
# 
# dat_sds <- dat %>%
#   na.omit() %>%
#   select(all_covs) %>%
#   summarise_all(sd)

# covariates_predict_clean <- covariates_predict %>%
#   as.data.frame() %>%
#   na.omit() %>%
#   mutate_at(vars(all_covs), scale(center = colMeans(cols), scale = apply(na.omit(cols), 2, sd))) %>%
#   mutate_at(vars(all_covs), as.numeric)

scale_covs <- function (covs, means, sds) {
  cols_scale <- match(names(means), colnames(covs))
  covs_sub <- covs[, cols_scale]
  covs_sub <- sweep(covs_sub, 2, means, "-")
  covs_sub <- sweep(covs_sub, 2, sds, "/")
  covs[, cols_scale] <- covs_sub
  covs
}

cols <- select(dat_clean, all_covs)
dat_means <- colMeans(cols)
dat_sds <- apply(cols, 2, sd)
covariates_predict <- scale_covs(raster::extract(covariates_cropped, seq(1, ncell(covariates_cropped), 1), df = TRUE), dat_means, dat_sds) %>%
  na.omit()

#prediction arrays
survival_predict <- get_survival(covariates_predict)
fecundity_predict <- get_fecundity(covariates_predict)
lambdas_predict <- get_lambda(survival_predict, fecundity_predict)
rate_predict <- get_rate(lambdas_predict)

# map posterior mean estimates
# idx is index of nonNA cells in prediction array
idx <- which(!is.na(getValues(covariates_cropped[[1]])))
  
map_variable <- function (greta_array, draws) {
  vals <- calculate(greta_array, draws)
  vals_mat <- as.matrix(vals)
  mean <- colMeans(vals_mat)
  map <- covariates_cropped[[1]]
  map[idx] <- mean
  map
}

lambda_map <- map_variable(lambdas_predict, draws)
rate_map <- map_variable(rate_predict, draws)
fec_map <- map_variable(fecundity_predict, draws)
surv_map <- map_variable(survival_predict$adult, draws)
