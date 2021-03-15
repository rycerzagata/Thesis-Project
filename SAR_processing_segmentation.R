library(raster)
library(rgrass7)
library(sp)
library(terra)
library(rgdal)

use_sp()

setwd("C:/TMP")
# setwd("C:/Users/HP/Documents/Thesis/Data")

# Set on-disk raster variable
rname <- paste(getwd(), "sumatra.tif", sep="/")

# Set GRASS environment and database location
loc <- initGRASS("C:/OSGEO4~1/apps/grass/grass78",
                 home=getwd(), gisDbase="GRASS_TEMP", override=TRUE )

# Import raster to GRASS and set region
execGRASS("r.in.gdal", flags="o", parameters=list(input=rname, output="tmprast"))
execGRASS("g.region", parameters=list(raster="tmprast") ) 

# Calculate 9x9 focal mean 
execGRASS("r.neighbors", flags="overwrite", parameters=list(input="tmprast", output="xxfm", 
                                                            method="average", size=as.integer(9)) )

img1 <- readRAST("tmprast")

# Computing annual temporal statistics of SAR data

img1 <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/MScAgata_S1_VH_TileID_1291.tif")
crs(img1) <- CRS("+init=epsg:31981")

# Compute temporal statistics of the SAR data for VH
img1_avg <- calc(img1, fun = function(x) {mean(x, na.rm=TRUE)})
img1_md <- calc(img1, fun = function(x) {median(x, na.rm=TRUE)})
img1_sd <- calc(img1, fun = function(x) {sd(x, na.rm=TRUE)})
img1_q <- calc(img1, fun = function(x) {quantile(x,probs = c(0.10, 0.90),na.rm=TRUE)} )
img1_q10 <- img1_q$layer.1
img1_q90 <- img1_q$layer.2
img1_stack <- stack(img1_avg, img1_md, img1_sd, img1_q10, img1_q90)

# Compute local statistics of the SAR data for VH
img1_avg_f <- focal(img1_avg, w=matrix(1/25,nrow=5,ncol=5), fun=mean)
img1_sd_f <- focal(img1_sd, w=matrix(1/25,nrow=5,ncol=5), fun=sd)

# Create a stack from all the statistics for VH image, repeat for VV or load the already calculated stats (see below)
SAR_stack <- stack(img1_avg, img1_md, img1_sd, img1_q10, img1_q90, img1_avg_f, img1_sd_f)
names(SAR_stack) <- c('avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh') # check the importance of the temporal/spatial features
writeRaster(SAR_stack, filename = "C:/Users/HP/Documents/Thesis/Data/Radar_test/cropped_SAR_stack_tile1291.tif", format= 'GTiff')

# Crop to an extent
SAR_stack_VV <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/SAR_stack_tile1291_VV.tif")
SAR_stack_VH <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/SAR_stack_tile1291_VH.tif")
r <- addLayer(SAR_stack_VV, SAR_stack_VH)
names(r) <- c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh')
e <- extent(c(735000, 745000, 9605000, 9615000))
plot(extent(r))
r <- crop(r, e)

# Spatial dimensions, local statistics
# 2 sets of features, one for VV and one for VH

#### Image segmentation (mean of the SAR image as an input) ####

# Smoothen the image
# Apply a gaussian blur with sigma=3
gf <- focalWeight(r, 3, "Gauss")
s1 <- focal(r$avg_vv, w=gf, na.rm=TRUE)
s2<- focal(r$md_vv, w=gf, na.rm=TRUE)
s3 <- focal(r$sd_vv, w=gf, na.rm=TRUE)
s4 <- focal(r$q10_vv, w=gf, na.rm=TRUE)
s5 <- focal(r$q90_vv, w=gf, na.rm=TRUE)
s6 <- focal(r$q90_vv, w=gf, na.rm=TRUE)
s7 <- focal(r$avg_vh, w=gf, na.rm=TRUE)
s8<- focal(r$md_vh, w=gf, na.rm=TRUE)
s9 <- focal(r$sd_vh, w=gf, na.rm=TRUE)
s10 <- focal(r$q10_vh, w=gf, na.rm=TRUE)
s11 <- focal(r$q90_vh, w=gf, na.rm=TRUE)
s12 <- focal(r$q90_vh, w=gf, na.rm=TRUE)
smoothen <- stack(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12)

writeRaster(smoothen, filename = "C:/Users/HP/Documents/Thesis/Data/Radar_test/cropped_SAR1291_stack_gauss_s12.tif", format= 'GTiff', overwrite=TRUE)

# Orfeo Toolbox segmentation - Large Scale Mean-Shift algorithm

system('cd C:/Users/HP/Documents/Thesis/Data/Radar_test')
system('"C:/Users/HP/Documents/Thesis/OTB/bin/otbcli_LargeScaleMeanShift" -in "C:/Users/HP/Documents/Thesis/Data/Radar_test/cropped_SAR1291_stack_gauss_s12.tif" -spatialr 25 -ranger 1.3 -minsize 20 -mode vector -mode.vector.out "C:/Users/HP/Documents/Thesis/Data/Radar_test/meanshift16.shp"')

# Read the output shp file into R
# MEANSHIFT 12 WAS THE BEST SO FAR
out_path = 'C:/Users/HP/Documents/Thesis/Data/Radar_test'
lyr_name = 'meanshift16'

# Check the segmentation results
shp <- readOGR(dsn=out_path, layer = lyr_name)
plot(shp)

# Orfeo Toolbox segmentation -  watershed algorithm

# system('C:/Users/HP/Documents/Thesis/OTB/bin/otbcli_Segmentation -in C:/Users/HP/Documents/Thesis/Data/Radar_test/cropped_SAR1291_stack_gauss_s3.tif -mode vector -mode.vector.out C:/Users/HP/Documents/Thesis/Data/Radar_test/watershed12.shp -filter watershed -filter.watershed.threshold 0.11 -filter.watershed.level 0.16')


#### Train and validate the random forest model ####

# Load the shapefile with segments overlapping with the ALS strip. The segments are labeled as forest/non-forest.
# The ALS and polygons were processed in ArcGIS with the output file called meanshift_16_labeled.shp.
out_path = 'C:/Users/HP/Documents/Thesis/Data'
lyr_name = 'meanshift16_labeled'
labeled  <- readOGR(dsn=out_path, layer = lyr_name)

# Get the coordinates of the centers of polygons.
library(rgeos)
centr <- gCentroid(labeled, byid = TRUE)
xy <- coordinates(centr)

# Prepare training/validation data for segments (polygon vector) that are overlapping the ALS strip
# Add the class label to the centroids table according the centroid/segment ID

library(spatialEco)
library(exactextractr)

v <- zonal.stats(labeled, r, stats = "mean") # compute stats for each polygon
v <- cbind (v, as.data.frame(xy))
v <- cbind(v, as.data.frame(labeled$gridcode)) # add labels to the df
names(v) <- c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh', 'x', 'y', 'class')
spLabeled <- v


# Prepare a random sample split between training and validation (70/30)
NumberOfRows <- nrow(spLabeled)
NumToExcludeForValidation <- (NumberOfRows%/%3.3333)
NumberForTraining <- (NumberOfRows - NumToExcludeForValidation)

index <- 1:NumberOfRows
ValidIndex <- sample(index, NumToExcludeForValidation)
trainingAll <- spLabeled[-ValidIndex,]
trainingAll <- na.omit(trainingAll)
training <- trainingAll[1:16]
trainingClass <- trainingAll[17]
trainingClass <- trainingClass[,1]
testAll <- spLabeled[ValidIndex,]
test <- testAll[1:16]
testClass <- testAll[17]
testClass <- testClass[,1]



# Random Forest
library(randomForest)
rfModel <- randomForest(x = training, y = trainingClass, na.action = na.omit, importance = TRUE)
rfPredict <- predict(rfModel, test)
library(caret)
confusionMatrix(rfPredict, testClass) # 97%

# Confusion Matrix and Statistics

#            Reference
# Prediction   0   1
#          0   8   1
#          1   8 190

# Accuracy : 0.9565          
# 95% CI : (0.9191, 0.9799)
# No Information Rate : 0.9227          
# P-Value [Acc > NIR] : 0.03757         

# Kappa : 0.6188    

###############################

# Use the new ML model to label the centroids (point vector) to either forest or non-forest

out_path = 'C:/Users/HP/Documents/Thesis/Data/Radar_test'
lyr_name = 'meanshift16'
polys  <- readOGR(dsn=out_path, layer = lyr_name)
# get centroids
library(rgeos)
centr <- gCentroid(polys, byid = TRUE)
xy <- coordinates(centr)


# Extract the values of SAR pixels within polygons
library(spatialEco)
library(exactextractr)
stats <- zonal.stats(polys, r, stats = "mean")
names(stats) <- c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh')
stats <- as.matrix(stats)
stats_coords <- cbind(stats, xy)
stats_coords <- na.omit(stats_coords)
coords <- stats_coords[, 15:16]
# Predict the class of segments based on the statistics extracted from SAR image within the polygons.
segmPredict <- predict(rfModel, stats_coords)
segmPredict <- cbind(stats_coords[,1:14], as.data.frame(segmPredict))
spSegmPred <- SpatialPointsDataFrame(coords , segmPredict, proj4string = CRS("+init=epsg:31981"))
writeOGR(spSegmPred, dsn=out_path, layer= 'segment_prediction_points', driver="ESRI Shapefile", overwrite_layer=TRUE)

# Transfer the centroid label to the corresponding segment (polygon vector) using spatial join in ArcGIS
# Colour the segments (polygon vector) according to the predicted label to visualize the results (rasterize the centroids using their predicted label to get the final raster forest/non-forest map).


 
 
 
