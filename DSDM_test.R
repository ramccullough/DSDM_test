#Install required packages
# install.packages("zoon")
# install.packages("devtools")
# install.packages("bayesplot")
# install.packages("ggplot2")

library(raster)
library(zoon)
library(greta)
library(bayesplot)
library(dplyr)
library(rgeos)
library(tidyr)
library(RColorBrewer)
library(ggplot2)

# 1.2 Set working directory for downloading data  -----------------------------
setwd("~/")
rel_path_input <- 'pp_grids/pp_grids/'

r <- aggregate(shapefile('pp_grids/pp_grids/cma.shp'), dissolve = TRUE) %>%
  as("SpatialPolygonsDataFrame")

# 2.4 Create raster stack of environmental covariates -------------------------
ndvi <- mask(raster(x = paste0(rel_path_input, 'ndvi.tif')), r) #wetness index
temp <- mask(raster(x = paste0(rel_path_input, 'meananntemp.tif')), r)
prec <- mask(raster(x = paste0(rel_path_input, 'annprecip.tif')), r)

covariates <- stack(temp,
                    prec,
                    ndvi)

missing <- calc(covariates, anyNA)
mask <- covariates[[1]]*0
mask[missing] <- NA
covariates <- mask(covariates, mask)

# crop
e <- new("Extent", xmin = 143.9, xmax = 145.8,
         ymin = -38.8, ymax = -37.7)
covariates_cropped <- crop(covariates, e)
covariates_cropped_matrix <- as.matrix(covariates_cropped)

# scaling values
covariate_means <- colMeans(covariates_cropped_matrix, na.rm = TRUE)
covariate_sds <- apply(covariates_cropped_matrix, 2, sd, na.rm = TRUE)

covariates_cropped <- scale(covariates_cropped, covariate_means, covariate_sds)

# Define covariates for fecundity and survival regression models
survival_covs <- c('meananntemp', 'annprecip', 'ndvi')
fecundity_covs <- c('meananntemp', 'annprecip')

all_covs <- c(survival_covs, fecundity_covs)

ncov_survival <- length(survival_covs)
ncov_fecundity <- length(fecundity_covs)

# extract environmental covariates for Melbourne bay locations
cells <- which(!is.na(getValues(covariates_cropped[[1]])))
n <- 500

sites <- sampleRandom(covariates_cropped, n, na.rm = TRUE, cells = TRUE, xy= TRUE)

vals <- raster::extract(covariates_cropped, sites[, "cell"])

vals <- as.matrix(vals)
dat <- vals[, all_covs]

# Survival and fecundity data for ringtail possums from McCarthy, Lindenmayer and Possingham
# adults survival =0.64
# juvenile survival = 0.33
# newborn females/breeder = 1.23
default_logit_survival_adult <- qlogis(0.64)
default_logit_survival_juvenile <- qlogis(0.33)
default_log_fecundity <- log(1.23)

# choose arbritrary regression coefficients for create simulated data

b_fecundity <- matrix(c(0.05555799, 0.02175256), 2, 1)
b_survival <- matrix(c(0.005308345, 0.002339500, -0.108931109), 3, 1)
l_intercept <- 0.1

## Priors for beta parameters
beta_fecundity <- normal(0, 0.1, dim = ncov_fecundity)
beta_survival <- normal(0, 0.1, dim = ncov_survival)
likelihood_intercept <- variable()

#survival
get_survival <- function(covs){
  x_survival <- dat[, covs]
  survival_logit_adult <- default_logit_survival_adult + x_survival %*% beta_survival
  survival_logit_juvenile <- default_logit_survival_juvenile + x_survival %*% beta_survival
  adult_survival <- ilogit(survival_logit_adult)
  juvenile_survival <- ilogit(survival_logit_juvenile)
  survival <- list(adult = adult_survival, juvenile = juvenile_survival)
}

#fecundity
get_fecundity <- function(covs){
  x_fecundity <- dat[, covs]
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

survival <- get_survival(survival_covs)
fecundity <- get_fecundity(fecundity_covs)
lambda <- get_lambda(survival, fecundity)

# poisson rate for abundance
pi <- exp(likelihood_intercept)*lambda
pi_true <- calculate(pi, values = list(likelihood_intercept = l_intercept, beta_fecundity = b_fecundity, beta_survival = b_survival))

# simulated abundance data
y <- rpois(n, pi_true)

# abundance distribution
distribution(y) <- poisson(pi)
  
m <- model(beta_fecundity, beta_survival, likelihood_intercept)

chains = 4
nsamples = 1000

draws <- greta::mcmc(m, n_samples = nsamples, chains = chains)

# summary stats
outputs <- list(
  r_hat = coda:::gelman.diag(draws, multivariate = FALSE),
  n_eff = coda:::effectiveSize(draws),
  summary = summary(draws))

# change mcmc.list variable names for plotting
varnames(draws) <- c("beta^{F[1]}", "beta^{F[2]}", "beta^{S[1]}", "beta^{S[2]}", "beta^{S[3]}", "alpha")


mcmc_data <- mcmc_trace_data(draws)

Blues <- c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")

# Trace plot

# Set folder for all output plots
setwd("~/Plots/Simulation")

png("trace.png",
    height = 800,
    width = 1200)
ggplot(data = mcmc_data, aes(x = iteration, y = value)) +
  geom_line(aes(colour = chain)) +
  facet_wrap(~parameter, labeller = label_parsed) +
  scale_colour_manual(name = "", values = c(Blues[[6]], Blues[[4]], Blues[[3]], Blues[[1]])) +
  guides(colour = guide_legend(title = "Chain")) +
  labs(x = "Iteration", y = "Parameter value") + 
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30, face="bold"),
        strip.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))
dev.off()
  
# Histogram
data_true <- data.frame(parameter = varnames(draws), true_value = c(0.05555799, 0.02175256, 0.005308345, 0.002339500, -0.108931109, 0.1))

png("histogram.png",
    height = 800,
    width = 1200)
ggplot() + 
  geom_histogram(data = mcmc_data, aes(x = value, y = ..density.., fill = "HMC draws"), colour = Blues[[3]]) +
  scale_fill_manual(name = "", values = Blues[[2]]) +
  facet_wrap(~parameter, labeller = label_parsed) +
  geom_vline(data = data_true, aes(xintercept = true_value, colour = "True value"), size = 1.0) +
  scale_colour_manual(name = "", values = Blues[[6]]) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  labs(x = "Parameter value", y = "Density") + 
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30, face="bold"),
        strip.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28),
        panel.spacing.x = unit(12, "mm"))
dev.off()

# PPC
pi_sim <- calculate(pi, draws)
pi_sim_matrix <- as.matrix(pi_sim)

y_sim <- array(NA, dim = dim(pi_sim_matrix))
y_sim[] <- rpois(length(pi_sim_matrix), pi_sim_matrix[])


ppc_data <- ppc_data(y, y_sim)

png("ppc_overlay.png",
    height = 800,
    width = 1000)
ggplot(ppc_data) +
  stat_ecdf(data = function(x) dplyr::filter(x, !.data$is_y), aes(x = value,  group = rep_id, colour = "Model predictions"), size = 0.25) +
  stat_ecdf(data = function(x) dplyr::filter(x, .data$is_y), aes(x = value, colour = "Observed data"), size = 1.0) +
  scale_colour_manual(name = "", values = c(Blues[[2]], Blues[[6]])) +
  guides(colour = guide_legend(title = "")) +
  labs(x = "Abundance value", y = "Cumulative density") +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30, face="bold"),
        strip.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))
dev.off()
  

# p-values
y_sim_jitter <- jitter(y_sim)
y_jitter <- jitter(y)
p <- colMeans(sweep(y_sim_jitter, 2, STATS = y_jitter, FUN = `>`))

p_data <- data.frame(quantile = p)

bin_breaks <- seq(0, 1, 0.1)

png("pvalue.png",
    height = 800,
    width = 1200)
ggplot(p_data) +
  geom_histogram(breaks = bin_breaks, aes(x = quantile), colour = Blues[[3]], fill = Blues[[2]]) +
  labs(x = "Posterior predictive p-value", y = "Count")  +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30, face="bold"),
        strip.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))
dev.off()


# QQ plot of quantile residuals
z_data <- data.frame(z_score = qnorm(p))

png("qqplot.png",
    height = 800,
    width = 800)
ggplot(z_data) +
  geom_qq(aes(sample = z_score), colour = Blues[[5]], size = 3) +
  geom_abline() +
  labs(x = "Standard normal quantiles", y = "Simulated data quantile residuals")  +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26, face="bold"),
        strip.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 24))
dev.off()

# plot simulated abundance data histogram
y_df <- data.frame(abundance = y)

png("example_sim_hist.png",
    height = 800,
    width = 1200)
ggplot(data = y_df) +
  geom_histogram(aes(x = abundance), binwidth = 1, colour = Blues[[3]], fill = Blues[[2]]) + 
  labs(x = "Simulated abundance (no. of individuals)", y = "Count")  +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30, face="bold"),
        strip.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 28))
dev.off()

# plot sites map with dots scaled by number present
sites_df <- as.data.frame(cbind(sites, y_df))
r_cropped <- crop(r, e)
map_cropped <- covariates_cropped[[1]]*0

plot(map_cropped,
     legend = FALSE,
     col = "white",
     axes = FALSE,
     box = FALSE)

plot(r_cropped, add = TRUE)

points(sites_df[sites_df[, 'abundance'] == 0, c('x', 'y')],
       pch = 1,
       cex = 0.8,
       col = Blues[[6]])

points(sites_df[sites_df[, 'abundance'] == 1, c('x', 'y')],
       pch = 20,
       cex = 0.8,
       col = Blues[[2]])

points(sites_df[sites_df[, 'abundance'] == 2, c('x', 'y')],
       pch = 20,
       cex = 1.2,
       col = Blues[[3]])

points(sites_df[sites_df[, 'abundance'] == 3, c('x', 'y')],
       pch = 20,
       cex = 1.6,
       col = Blues[[4]])

points(sites_df[sites_df[, 'abundance'] == 4, c('x', 'y')],
       pch = 20,
       cex = 2,
       col = Blues[[5]])

points(sites_df[sites_df[, 'abundance'] == 5, c('x', 'y')],
       pch = 20,
       cex = 2.4,
       col = Blues[[6]])

legend("right", legend = c("0", "1", "2", "3", "4", "5"), 
       bty = "n",
       col = c(Blues[[6]], Blues[[2]], Blues[[3]], Blues[[4]], Blues[[5]], Blues[[6]]), 
       pch = c(1, 20, 20, 20, 20, 20),
       pt.cex = c(0.8, 0.8, 1.2, 1.6, 2, 2.4))

