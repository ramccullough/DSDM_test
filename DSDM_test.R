#Install required packages
install.packages("zoon")
install.packages("devtools")
install.packages("bayesplot")

library(raster)
library(zoon)
library(greta)
library(bayesplot)
library(dplyr)

#set working directory for data
setwd("C:/Users/racha/Google Drive (rmccullough@student.unimelb.edu.au)/MSc/Research/Computational")

#Define functions
#scale covariates
scale_covs <- function (covs, means, sds) {
  cols_scale <- match(names(means), colnames(covs))
  covs_sub <- covs[, cols_scale]
  covs_sub <- sweep(covs_sub, 2, means, "-")
  covs_sub <- sweep(covs_sub, 2, sds, "/")
  covs[, cols_scale] <- covs_sub
  covs
}

LoadModule('Bioclim')
resolution <- 2.5

bioclim <- getData('worldclim', var = 'bio', res = resolution)

#Get bio data for Melbourne region
e <- new("Extent", xmin = 143.612155685872, xmax = 145.549489975136, 
         ymin = -38.8184069791436, ymax = -37.481080244538)
bioclim_cropped <- crop(bioclim, e)
#plot(bioclim_cropped[[1]], maxpixels = 100000)

#choose survival and fecundity covariates (to refine/edit later with ref to biology)
#survival
#bio1
#bio12
#bio6

#fecundity
#bio1
#bio12

survival_covs<-c('bio1', 'bio12', 'bio6')
fecundity_covs<-c('bio1', 'bio12')
all_covs<-c(survival_covs, fecundity_covs)

ncov_survival<-length(survival_covs)
ncov_fecundity<-length(fecundity_covs)

#extract environmental covariates for Melbourne bay locations
cells <- which(!is.na(getValues(bioclim_cropped[[1]])))
n<- length(cells)

sites<-sampleRandom(bioclim_cropped, size=n, na.rm=TRUE, cells=TRUE)

vals<-extract(bioclim_cropped, sites[, "cell"])

vals<- as.matrix(vals)
dat<-vals[, all_covs]

dat_means <- colMeans(dat)
dat_sds <- apply(dat, 2, sd)

dat<-scale_covs(dat, dat_means, dat_sds)

#survival and fecundity data for ringtail possums from McCarthy, Lindenmayer and Possingham
#adults survival =0.64
#juvenile survival =0.33
#newborn females/breeder =1.23
default_logit_survival_adult<-qlogis(0.64)
default_logit_survival_juvenile <- qlogis(0.33)
default_log_fecundity<-log(1.23)

#choose arbritrary regression coefficients for create simulated data
# b_fecundity<-as.matrix(rnorm(ncov_fecundity, mean=0, sd=0.1))
# b_survival<-as.matrix(rnorm(ncov_survival, mean=0, sd=0.1))
b_fecundity<-matrix(c(0.018, -0.095), 2, 1)
b_survival<-matrix(c(-0.183, -0.010, 0.036), 3, 1)
l_intercept<-0.1

##Using greta arrays
beta_fecundity <- normal(0, 1/ncov_fecundity, dim=ncov_fecundity)
beta_survival <- normal(0, 1/ncov_survival, dim=ncov_survival)
likelihood_intercept <- variable()

#survival
get_survival<-function(covs){
  x_survival<-dat[, covs]
  survival_logit_adult<-default_logit_survival_adult+x_survival%*%beta_survival
  survival_logit_juvenile<-default_logit_survival_juvenile+x_survival%*%beta_survival
  adult_survival<-ilogit(survival_logit_adult)
  juvenile_survival<-ilogit(survival_logit_juvenile)
  survival<-list(adult=adult_survival, juvenile=juvenile_survival)
}

#fecundity
get_fecundity<-function(covs){
  x_fecundity<-dat[, covs]
  fecundity_log=default_log_fecundity+x_fecundity%*%beta_fecundity
  fecundity=exp(fecundity_log)
}

#lambda
get_lambda<-function(survival, fecundity){
  top_row <- cbind(fecundity * survival$juvenile,
                   fecundity * survival$adult)
  bottom_row <- cbind(survival$juvenile,
                      survival$adult)
  matrices <- abind(top_row, bottom_row, along = 3)
  iterated <- greta.dynamics::iterate_matrix(matrices, niter = 20)
  lambda<-iterated$lambda
}

survival<-get_survival(survival_covs)
fecundity<-get_fecundity(fecundity_covs)
lambda<-get_lambda(survival, fecundity)

#poisson rate for abundance
# pi<-exp(likelihood_intercept)*lambda
# pi_true<-calculate(pi, values=list(likelihood_intercept=l_intercept, beta_fecundity=b_fecundity, beta_survival=b_survival))

#probability of presence
p<-icloglog(likelihood_intercept+log(lambda))
p_true<-calculate(p, values=list(likelihood_intercept=l_intercept, beta_fecundity=b_fecundity, beta_survival=b_survival))

#simulated PA data
y<-rbinom(n, 1, p_true)

#simulated abundance data
# y<-rpois(n, pi_true)

#PA distribution
distribution(y)<-bernoulli(p)

#abundance distribution
# distribution(y)<-poisson(pi)

m<-model(beta_fecundity, beta_survival, likelihood_intercept)
#plot(m)

chains=4
nsamples=1000

draws<-mcmc(m, n_samples=nsamples, chains = chains)
summary(draws)
mcmc_hist(draws)
mcmc_recover_hist(draws, true=c(b_fecundity, b_survival, l_intercept))


# o<-opt(m, hessian=TRUE)
#  
# fec_draws <- calculate(beta_fecundity, draws)
# fec_matrix <- as.matrix(draws)
# cor(fec_matrix)
# 
# surv_draws <- calculate(beta_survival, draws)
# surv_matrix <- as.matrix(draws)
# cor(surv_matrix)
#  
# sqrt(solve(h$hessian$beta_fecundity))

#get mle and covariances matrices
summary(glm(y~lambda, family=stats::binomial()))

#find negative log likelihood
nll<-function(){
  p<-icloglog(likelihood_intercept+log(lambda))
  -sum(y*log(p)+(1-y)*log(1-p))
}



