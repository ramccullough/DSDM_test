# =============================================================================
# Daily maximum temperature data for Victoria sourced from 
# https://data.csiro.au/dap/landingpage?execution=e1s2
#
# Aim of script below is to load, clean and extract relevant data. Final output
# is to be a raster of the average number of days per year where the maximum 
# temperature is above 36C for the Port Phillip Bay area. Average is to be 
# taken from 10 most recent years of data (2004 - 2014).
# =============================================================================

install.packages('ncdf4')
install.packages('raster')
library(ncdf4)
library(raster)

# 1. Open file as raster ------------------------------------------------------
file <- '/Users/rmccullough/Documents/Research/daily_maxt_3ds_19810101_20140731.nc'
# file_nc <- nc_open(file)
# file_raster <- stack(file_nc)

file_raster <- brick(file)

e <- new("Extent", xmin = 143.9, xmax = 145.8, 
         ymin = -38.8, ymax = -37.7)

# 2. Extract layers for X2004.07.31 to X2014.07.31 (10 years of data) ---------
raster_10y <- file_raster[[8613:12265]]

# 3. Extract one year (365 layers) of values for testing ----------------------
raster_y1 <- raster_10y[[1:365]]
raster_y1_cropped <- crop(raster_10y, e) # not working/taking too long

# 3a. Extract other years of values -------------------------------------------
raster_y2 <- raster_10y[[366:730]]
raster_y3 <- raster_10y[[731:1095]]
raster_y4 <- raster_10y[[1096:1461]]
raster_y5 <- raster_10y[[1462:1826]]
raster_y6 <- raster_10y[[1827:2191]]
raster_y7 <- raster_10y[[2192:2556]]
raster_y8 <- raster_10y[[2557:2922]]
raster_y9 <- raster_10y[[2923:3287]]
raster_y10 <- raster_10y[[3288:3652]]
  
# 4. Extract one day (1 layer) of values for testing --------------------------
raster_1d <- raster_y1[[1]]

# 4. Aggregrate to reduce size for testing ------------------------------------
raster_10y <- aggregate(raster_10y, 10) # currently not working/taking too long
raster_y1 <- aggregate(raster_y1, 10) # currently not working/taking too long
raster_1d <- aggregate(raster_1d, 10) # works

# 5. Test plot of single layer ------------------------------------------------
raster_1d_cropped <- crop(raster_1d, e)
plot(raster_test_cropped)

# 6. Find and mask missing values ---------------------------------------------
calc(raster_y1, anyNA) # currently returning an error

# Testing with 1 later
calc(raster_1d, anyNA) # works 

# Testing with 2 layers
raster_2d <- raster_y1[[1:2]]
missing <- calc(raster_2d, anyNA) # returns error 
# mask <- raster_2d[[1]]*0
# mask[missing] <- NA
# raster_2d <- mask(raster_2d, mask)

# 7. Save 10 years of cropped data as raster ----------------------------------
setwd('/Users/rmccullough/Documents/Research')
# writeRaster(raster_extract, 'vic_daily_max_2004_2014_raster.tif', format = 'GTiff) 
# not working/taking too long

# Save one year of data
# writeRaster(raster_1y, 'vic_daily_max_1y.tif', format = "GTiff")
# not working/taking too long

# Saving as separate layers (Tiff format)
writeRaster(raster_y1, 'vic_daily_max_1y/layer.tif', format = "GTiff", bylayer = TRUE)

#7a Saving other years of data as separate layers -----------------------------
writeRaster(raster_y2, 'vic_daily_max_y2/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y3, 'vic_daily_max_y3/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y4, 'vic_daily_max_y4/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y5, 'vic_daily_max_y5/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y6, 'vic_daily_max_y6/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y7, 'vic_daily_max_y7/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y8, 'vic_daily_max_y8/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y9, 'vic_daily_max_y9/layer.tif', format = "GTiff", bylayer = TRUE)
writeRaster(raster_y10, 'vic_daily_max_y10/layer.tif', format = "GTiff", bylayer = TRUE)

# 8. Re-open saved TIF files  -------------------------------------------------
# Test 10 days of data
# How to simplify this for more layers?
# d1 <- raster('vic_daily_max_1y_1.tif')
# d2 <- raster('vic_daily_max_1y_2.tif')
# d3 <- raster('vic_daily_max_1y_3.tif')
# d4 <- raster('vic_daily_max_1y_4.tif')
# d5 <- raster('vic_daily_max_1y_5.tif')
# d6 <- raster('vic_daily_max_1y_6.tif')
# d7 <- raster('vic_daily_max_1y_7.tif')
# d8 <- raster('vic_daily_max_1y_8.tif')
# d9 <- raster('vic_daily_max_1y_9.tif')
# d10 <- raster('vic_daily_max_1y_10.tif')

# Stack test layers
# How to simplify this for more layers?
# raster_10d <- stack(d1,
#                     d2,
#                     d3,
#                     d4,
#                     d5,
#                     d6,
#                     d7,
#                     d8,
#                     d9,
#                     d10)
# 
# missing <- calc(raster_10d, anyNA)
# mask <- raster_10d[[1]]*0
# mask[missing] <- NA
# raster_10d <- mask(raster_10d, mask)

# Need to repeat step 8 for remaining 9 years of data 
# Layer names
layer_names_y1 <- sapply('vic_daily_max_y1/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y2 <- sapply('vic_daily_max_y2/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y3 <- sapply('vic_daily_max_y3/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y4 <- sapply('vic_daily_max_y4/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y5 <- sapply('vic_daily_max_y5/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y6 <- sapply('vic_daily_max_y6/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y7 <- sapply('vic_daily_max_y7/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y8 <- sapply('vic_daily_max_y8/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y9 <- sapply('vic_daily_max_y9/layer_', paste, 1:365, ".tif", sep = "")
layer_names_y10 <- sapply('vic_daily_max_y10/layer_', paste, 1:365, ".tif", sep = "")

# Stack layers for each year
raster_y1 <- stack(sapply(layer_names_y1, raster))
raster_y2 <- stack(sapply(layer_names_y2, raster))
raster_y3 <- stack(sapply(layer_names_y3, raster))
raster_y4 <- stack(sapply(layer_names_y4, raster))
raster_y5 <- stack(sapply(layer_names_y5, raster))
raster_y6 <- stack(sapply(layer_names_y6, raster))
raster_y7 <- stack(sapply(layer_names_y7, raster))
raster_y8 <- stack(sapply(layer_names_y8, raster))
raster_y9 <- stack(sapply(layer_names_y9, raster))
raster_y10 <- stack(sapply(layer_names_y10, raster))

# Define masking function
mask_missing <- function(stack){
  missing <- calc(stack, anyNA)
  mask <- stack[[1]]*0
  mask[missing] <- NA
  stack <- mask(stack, mask)
  stack_cropped <- crop(stack, e)
  return(stack_cropped)
}

raster_clean_y1 <- mask_missing(raster_y1)
raster_clean_y2 <- mask_missing(raster_y2)
raster_clean_y3 <- mask_missing(raster_y3)
raster_clean_y4 <- mask_missing(raster_y4)
raster_clean_y5 <- mask_missing(raster_y5)
raster_clean_y6 <- mask_missing(raster_y6)
raster_clean_y7 <- mask_missing(raster_y7)
raster_clean_y8 <- mask_missing(raster_y8)
raster_clean_y9 <- mask_missing(raster_y9)
raster_clean_y10 <- mask_missing(raster_y10)

# currently not working/taking too long

# 9. For each year, add up how many days with max temp over 36C at each site --
# Test with 10 days of data
max_days <- sum(raster_10d > 36)

setwd('/Users/rmccullough/Documents/Research/yearly_count')
max_days_y1 <- sum(raster_y1 > 36)
writeRaster(max_days_y1, 'max_days_y1.tif', format = "GTiff")

max_days_y2 <- sum(raster_y2 > 36)
writeRaster(max_days_y2, 'max_days_y2.tif', format = "GTiff")

max_days_y3 <- sum(raster_y3 > 36)
writeRaster(max_days_y3, 'max_days_y3.tif', format = "GTiff")

max_days_y4 <- sum(raster_y4 > 36)
writeRaster(max_days_y4, 'max_days_y4.tif', format = "GTiff")

max_days_y5 <- sum(raster_y5 > 36)
writeRaster(max_days_y5, 'max_days_y5.tif', format = "GTiff")

max_days_y6 <- sum(raster_y6 > 36)
writeRaster(max_days_y6, 'max_days_y6.tif', format = "GTiff")

max_days_y7 <- sum(raster_y7 > 36)
writeRaster(max_days_y7, 'max_days_y7.tif', format = "GTiff")

max_days_y8 <- sum(raster_y8 > 36)
writeRaster(max_days_y8, 'max_days_y8.tif', format = "GTiff")

max_days_y9 <- sum(raster_y9 > 36)
writeRaster(max_days_y9, 'max_days_y9.tif', format = "GTiff")

max_days_y10 <- sum(raster_y10 > 36)
writeRaster(max_days_y10, 'max_days_y10.tif', format = "GTiff")

# 10. Find average number of days per year over 10 years ----------------------
max_days_y1 <- raster('max_days_y1.tif')
max_days_y2 <- raster('max_days_y2.tif')
max_days_y3 <- raster('max_days_y3.tif')
max_days_y4 <- raster('max_days_y4.tif')
max_days_y5 <- raster('max_days_y5.tif')
max_days_y6 <- raster('max_days_y6.tif')
max_days_y7 <- raster('max_days_y7.tif')
max_days_y8 <- raster('max_days_y8.tif')
max_days_y9 <- raster('max_days_y9.tif')
max_days_y10 <- raster('max_days_y10.tif')

max_days_10y <- stack(max_days_y1,
                      max_days_y2,
                      max_days_y3,
                      max_days_y4,
                      max_days_y5,
                      max_days_y6,
                      max_days_y7,
                      max_days_y8,
                      max_days_y9,
                      max_days_y10)

av_days_data_10y <- mean(max_days_10y)
writeRaster(av_days_data_10y , 'av_days_data_10y.tif', format = "GTiff")

# =============================================================================
# END OF SCRIPT
# =============================================================================

# Testing space ---------------------------------------------------------------
# # list of layer names
# layer_names <- sapply('vic_daily_max_y1/layer_', paste, 1:365, ".tif", sep = "")
# 
# raster_stack <- stack(sapply(layer_names, raster))
# 
# # missing <- calc(raster_stack, anyNA)
# # mask <- raster_stack[[1]]*0
# # mask[missing] <- NA
# # raster_stack <- mask(raster_stack, mask)
# 
# max_days <- sum(raster_stack > 36)
# 
# raster_test <- raster('')