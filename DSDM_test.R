install.packages("rlang")
install.packages("zoon")
install.packages("devtools")
install.packages("abind")
install.packages("VGAM")

library(raster)
library(zoon)
library(greta)
library(abind)
library(VGAM)
LoadModule('Bioclim')

bioclim<-getData('worldclim', var='bio', res=10)

#Get bio data for Melbourne region
e <- new("Extent", xmin = 143.612155685872, xmax = 145.549489975136, 
         ymin = -38.8184069791436, ymax = -37.481080244538)
bioclim_cropped <- crop(bioclim, e)
plot(bioclim_cropped[[1]], maxpixels = 100000)

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

#extract environmental covariates for Melbourne bay locations

#sample_sites<-sampleRandom(bioclim_cropped, size=50, cells=TRUE)
#n<- ncell(bioclim_cropped)
#sites<-seq(1, ncell(bioclim_cropped), 1)

sites <- which(!is.na(getValues(bioclim_cropped[[1]])))
n<- length(sites)
vals<-extract(bioclim_cropped, sites)

vals<- as.matrix(vals)
dat<-vals[, all_covs]

#scale covariates
scale_covs <- function (covs, means, sds) {
  cols_scale <- match(names(means), colnames(covs))
  covs_sub <- covs[, cols_scale]
  covs_sub <- sweep(covs_sub, 2, means, "-")
  covs_sub <- sweep(covs_sub, 2, sds, "/")
  covs[, cols_scale] <- covs_sub
  covs
}

dat_means <- colMeans(dat)
dat_sds <- apply(dat, 2, sd)

dat<-scale_covs(dat, dat_means, dat_sds)


#survival and fecundity data for ringtail possums from McCarthy, Lindenmayer and Possingham
#adults survival =0.64
#juvenile survival =0.33
#newborn females/breeder =1.23
#sd for demographic parameters??
# default_logit_survival_adult<-normal(qlogis(0.64), 0.01 )
# default_logit_survival_juvenile <- normal(qlogis(0.33), 0.01)
# default_log_fecundity<-normal(log(1.23), 0.25)

##Using matrices

default_logit_survival_adult<-as.matrix(rep(qlogis(0.64), n))
default_logit_survival_juvenile <- as.matrix(rep(qlogis(0.33), n))
default_log_fecundity<-as.matrix(rep(log(1.23), n))

#choose arbritrary regression coefficients for create simulated data
ncov_fecundity<-length(fecundity_covs)
ncov_survival<-length(survival_covs)

beta_fecundity<-as.matrix(rnorm(ncov_fecundity, mean=0, sd=0.1))
beta_survival<-as.matrix(rnorm(ncov_survival, mean=0, sd=0.1))
likelihood_intercept<-0.1

#survival model
x_survival<-dat[, survival_covs]

survival_logit_adult<-default_logit_survival_adult+x_survival%*%beta_survival
survival_logit_juvenile<-default_logit_survival_juvenile+x_survival%*%beta_survival
adult_survival<-plogis(survival_logit_adult)
juvenile_survival<-plogis(survival_logit_juvenile)
survival<-list(adult=adult_survival, juvenile=juvenile_survival)

#fecundity model
x_fecundity<-dat[, fecundity_covs]

fecundity_log=default_log_fecundity+x_fecundity%*%beta_fecundity
fecundity=exp(fecundity_log)

#get lambda (intrinsic growth) from Leslie matrices
top_row <- cbind(fecundity * survival$juvenile,
                 fecundity * survival$adult)
bottom_row <- cbind(survival$juvenile,
                    survival$adult)
matrices <- abind(top_row, bottom_row, along = 3)

#cell-by-cell eigenvalues
lambdas<-array(0, dim=n)
for (i in 1:70){
  top<-matrices[i, ,1]
  bottom<-matrices[i, ,2]
  leslie_matrix<-abind(top, bottom, along=2)
  eigenvalues<-eigen(leslie_matrix)
  lambda<-eigenvalues$values
  lambdas[i]<-max(lambda)
}


#probability of presence
p<-1-exp(-(likelihood_intercept*lambdas))
y<-rbinom(n, 1, p)


##Using greta arrays
# b_fecundity <- normal(0, 1, dim=ncov_fecundity)
# b_survival <- normal(0, 1, dim=ncov_survival)
# l_intercept <- variable()
# 
# #survival
# survival_logit_adult<-default_logit_survival_adult+x_survival%*%b_survival
# survival_logit_juvenile<-default_logit_survival_juvenile+x_survival%*%b_survival
# adult_survival<-ilogit(survival_logit_adult)
# juvenile_survival<-ilogit(survival_logit_juvenile)
# 
# adult_survival_true<-calculate(adult_survival, values=list(b_survival=beta_survival))
# juvenile_survival_true<-calculate(juvenile_survival, values=list(b_survival=beta_survival))
# survival=list(adult=adult_survival_true, juvenile=juvenile_survival_true)
# 
# #fecundity
# fecundity_log=default_log_fecundity+x_fecundity%*%b_fecundity
# fecundity=exp(fecundity_log)
# fecundity<-calculate(fecundity, values=list(b_fecundity=beta_fecundity))
# 
# #get lambda from leslie matrices
# top_row <- cbind(fecundity * survival$juvenile,
#                fecundity * survival$adult)
# bottom_row <- cbind(survival$juvenile,
#                   survival$adult)
# matrices <- abind(top_row, bottom_row, along = 3)
# iterated <- greta.dynamics::iterate_matrix(matrices, niter = 10)
# 
# lambda<-iterated$lambda

#probability of presence
#p<-icloglog(likelihood_intercept+log(lambda))

#y<-rbinom(n, 1, p)





