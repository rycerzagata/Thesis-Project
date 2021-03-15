# Loading the required libraries
library(sp)
library(rgl)
library(raster)
library(rgeos)
library(bfastSpatial) # if it's not possible to install the package follow the commented code in the main script



# Setting working directory
setwd('C:/Users/HP/Documents/Thesis/Data/Fragstats')

# Change the path /filename of your map
FNFMap <- raster('TanDEM2.tif')

## Make an NA-value raster based on the LC raster attributes
ForestMask <- setValues(raster(FNFMap), NA)
## Assign 1 to all cells corresponding to the forest class
ForestMask[FNFMap==1] <- 1

ForestClumpsSieve <- areaSieve(ForestMask,thresh = 5000, directions = 8)
ForestClumpsSieve[is.na(ForestClumpsSieve)] <- 0

op <- par(mfrow=c(1, 2))
plot(FNFMap, legend=FALSE)
plot(ForestClumpsSieve, legend=FALSE)
par(op)

# Chenge the output file name
writeRaster(ForestClumpsSieve, filename = "C:/Users/HP/Documents/Thesis/Data/Global FNF maps/TanDEM_no_clumps.tif", format= 'GTiff')


