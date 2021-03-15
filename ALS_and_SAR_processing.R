# Loading the required libraries
library(lidR)
library(sp)
library(rgl)
library(TreeLS)
library(rLiDAR)
library(raster)
library(rgeos)
library(caret)


#### CHM COMPUTATION ####
# Setting working directory
setwd('C:/Users/HP/Documents/Thesis/Data/LiDAR_test/Merged')

# Preprocessing using LASTools software (doesn't have to be rerun)
system('C:/Users/HP>set PATH =%PATH%; C:/LAStools/bin')
system('C:/Users/HP>set RAW_FLIGHT_LINES= C:/Users/HP/Documents/Thesis/Data/LiDAR')
system('C:/Users/HP>set RAW_FORMAT=las')
system('C:/Users/HP>set TEMP_FILES=C:/lastools_temp')
system('C:/Users/HP>set OUTPUT_FILES=C:/lastools_output')
system('C:/LAStools/bin/lasboundary 
       -i C:/Users/HP/Documents/Thesis/Data/LiDAR/*.las -merged -o C:/lastools_output/merged.shp')

system('C:/LAStools/bin/lasground -i C:/Users/HP/Documents/Thesis/Data/LiDAR_test/*.las')
system('C:/LAStools/bin/lasheight -i C:/Users/HP/Documents/Thesis/Data/LiDAR_test/*.las')
system('C:/LAStools/bin/lasclassify -i C:/Users/HP/Documents/Thesis/Data/LiDAR_test/Singles/*.las')

system('C:/LAStools/bin/lasinfo -i C:/Users/HP/Documents/Thesis/Data/LiDAR_test/*.las')

system('C:/LAStools/bin/lasmerge -i C:/Users/HP/Documents/Thesis/Data/LiDAR_test/Singles/*.las -o C:/Users/HP/Documents/Thesis/Data/LiDAR_test/merged/merged3.las')

# Load the merged point cloud
lasfile <- "merged3.las"
memory.limit(5000000)

# Compute Canopy Height Model (CHM)
las1 <- readTLS(lasfile)
DTM <- grid_terrain(las1, res = 1, algorithm = knnidw(k=6L, p=2), keep_lowest = FALSE)
DSM <- normalize_height(las1, DTM) # OLD NAME OF THE FUNCTION: lasnormalize()
epsg(DSM) <- 31981
CHM <- grid_canopy(DSM, res=1, p2r(0.2))
CHM <- focal(CHM,w=matrix(1/9, nc=3, nr=3), na.rm=TRUE)
plot(CHM)


# Classify areas with >30% tree cover and trees higher than 5 meters as a forest
CHM[CHM > 5] <- 1
CHM[CHM < 1] <- 0
CHM[CHM > 1] <- 0
hist(CHM)
crs(CHM) <- CRS("+init=epsg:31981")
forestMap <- aggregate(CHM, fact = 10, fun=sum)
forestMap[forestMap > 30] <- 1
forestMap[forestMap < 1] <- 0
forestMap[forestMap > 1] <- 0

# Now remove the single cells (forest clumps) from the forest map using functions by BFAST Spatial package (use the commented code to download the package)

#library(remotes)
#Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
#remotes::install_github('loicdtx/bfastSpatial')
# Load the package
library(bfastSpatial)

# Make an NA-value raster based on the LC raster attributes
formask <- setValues(raster(forestMap), NA)
# Assign 1 to all cells corresponding to the forest class
formask[forestMap==0] <- 1

clumps <- areaSieve(formask, thresh = 5000, directions = 8)
clumps[is.na(clumps)] <- 10
op <- par(mfrow=c(1, 2))
plot(forestMap, legend=FALSE)
plot(clumps, legend=FALSE)
par(op)

## Make an NA-value raster based on the LC raster attributes
nonformask <- setValues(raster(forestMap), NA)
## Assign 0 to all cells corresponding to the non-forest class
nonformask[clumps==0] <- 0

forestMap_noclumps <- areaSieve(nonformask,thresh = 5000, directions = 8, keepzeros = TRUE)
forestMap_noclumps[is.na(forestMap_noclumps)] <- 1
plot(forestMap_noclumps)

# writeRaster(forestMap_noclumps, filename = "C:/Users/HP/Documents/Thesis/Data/forest_nonforest_map_CHM.tif", format = 'GTiff', overwrite = TRUE)
forestMap <- raster("C:/Users/HP/Documents/Thesis/Data/forest_nonforest_map_CHM.tif")

# Resample forest map using SAR image
SARimg <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/MScAgata_S1_VV_TileID_1291.tif")
crs(SARimg) <- CRS('+init=epsg:31981')
forestMap_resampled <- resample(forestMap, SARimg, method = 'ngb')

# Rerun everything for VH&VV

# Compute temporal statistics of the SAR data for VH
SARimg_avg <- calc(SARimg, fun = function(x) {mean(x, na.rm=TRUE)})
SARimg_md <- calc(SARimg, fun = function(x) {median(x, na.rm=TRUE)})
SARimg_sd <- calc(SARimg, fun = function(x) {sd(x, na.rm=TRUE)})
SARimg_q <- calc(SARimg, fun = function(x) {quantile(x,probs = c(0.10, 0.90),na.rm=TRUE)} )
SARimg_q10 <- SARimg_q$layer.1
SARimg_q90 <- SARimg_q$layer.2

# Compute local statistics of the SAR data for VH
SARimg_avg_f <- focal(SARimg_avg, w=matrix(1/25,nrow=5,ncol=5), fun=mean, na.rm=TRUE)
SARimg_sd_f <- focal(SARimg_sd, w=matrix(1/25,nrow=5,ncol=5), fun=sd, na.rm=TRUE)


# Create a stack from all the statistics for VH image
SAR_stack <- stack(SARimg_avg, SARimg_md, SARimg_sd, SARimg_q10, SARimg_q90, SARimg_avg_f, SARimg_sd_f)
# check the importance of the temporal/spatial features
writeRaster(SAR_stack, filename = "C:/Users/HP/Documents/Thesis/Data/Radar_test/SAR_stack_tile1291_VV.tif", format= 'GTiff')

SAR_stack_VV <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/SAR_stack_tile1291_VV.tif")
SAR_stack_VH <- brick("C:/Users/HP/Documents/Thesis/Data/Radar_test/SAR_stack_tile1291_VH.tif")
SAR_stack <- addLayer(SAR_stack_VV, SAR_stack_VH)
names(SAR_stack) <- c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh')

# Choose the areas within the ALS extent 
covs <- mask(SAR_stack, forestMap_resampled)
names(forestMap_resampled) <- "class"
trainingBrick <-addLayer(covs, forestMap_resampled)

# Extract all values into a matrix
valuetable <- getValues(trainingBrick)
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)
head(valuetable, n = 10)
valuetable$class <- factor(valuetable$class, levels = c(0:1))


#### Classification ####

# Split the data into training and test datasets with the 70/30 ratio
NumberOfRows <- nrow(valuetable)
NumToExcludeForValidation <- (NumberOfRows%/%3.3333)
NumberForTraining <- (NumberOfRows - NumToExcludeForValidation)

index <- 1:nrow(valuetable)
ValidIndex <- sample(index, NumToExcludeForValidation)
trainingAll <- valuetable[-ValidIndex,]
training <- trainingAll[1:14]
trainingClass <- trainingAll[15]
trainingClass <- trainingClass[,1]
testAll <- valuetable[ValidIndex,]
test <- testAll[1:14]
testClass <- testAll[15]
testClass <- testClass[,1]

## Test different machine learning classifiers

# Logistic Regression
library(stats)
glmModel <- glm(class ~ ., family = binomial(link = "logit"), data = trainingAll)
glmPredict <- predict(SAR_stack, model=glmModel, na.rm=TRUE, type="response") # type = "response" to get the probability
# the predict function gives the log-odds for a binomial function as output. You can also specify that you want to get the probability 
# (which is maybe easier to interprete). You can then maybe use a 50% probability threshold to distinguish between 0 and 1.
glmPredict[glmPredict >= 0.5] <- 1
glmPredict[glmPredict < 0.5] <- 0
writeRaster(glmPredict, filename = "C:/Users/HP/Documents/Thesis/Data/gmlForestMap.tif", format= 'GTiff', overwrite = TRUE)

# K-Nearest Neighbour
# Disadvantage: too slow when more variables are added
library(caret)
knnModel <- train(training, trainingClass, method = "knn")
knnPredict <- predict(SAR_stack, model=knnModel, na.rm=TRUE)
writeRaster(knnPredict, filename = "C:/Users/HP/Documents/Thesis/Data/knnForestMap.tif", format= 'GTiff', overwrite = TRUE)

# Random Forest
library(randomForest)
rfModel <- randomForest(class ~ ., data = trainingAll, importance = TRUE)
rfPredict <- predict(SAR_stack, model=rfModel, na.rm=TRUE)
varImpPlot(rfModel)
writeRaster(rfPredict, filename = "C:/Users/HP/Documents/Thesis/Data/rfForestMap.tif", format= 'GTiff', overwrite=TRUE)

library(e1071)
svmModel <- svm(class ~ ., data = trainingAll)
svmPredict <- predict(SAR_stack, model=svmModel, na.rm=TRUE)
writeRaster(svmPredict, filename = "C:/Users/HP/Documents/Thesis/Data/svmForestMap.tif", format= 'GTiff')

#### Validation of the ML models ####

# Create a confusion matrix for each method

pr1 <- predict(glmModel, test, type="response")
pr1[pr1 >= 0.5] <- 1
pr1[pr1 < 0.5] <- 0
prf <- factor(pr1, levels = c(0,1))
confusionMatrix(prf, testClass) # Accuracy : 0.8562 # Accuracy for VV and VH: 0.923

pr2 <- predict(knnModel, test)
confusionMatrix(pr2, testClass) # Accuracy : 0.9073 # Accuracy for VV and VH: 0.9157 

pr3 <- predict(rfModel, test)
confusionMatrix(pr3, testClass) #   Accuracy for VV only : 0.9146 Accuracy for VV and VH: 0.9397

pr4 <- predict(svmModel, test)
confusionMatrix(pr4, testClass) # Accuracy : 0.9112 Accuracy for VV and VH: 0.9251
# at the end come back to this and check the values using test dataset with x and y columns


rfForestMap <- raster("C:/Users/HP/Documents/Thesis/Data/rfForestMap.tif")

## Choose the best prediction based on calculated overall accuracy
## Now remove single cells (clumps) from the map using tools from BFAST Spatial package

# Make an NA-value raster based on the LC raster attributes
formask <- setValues(raster(rfForestMap), NA)
# Assign 1 to all cells corresponding to the forest class
formask[rfForestMap==1] <- 1

forclumpsSieve <- areaSieve(formask,thresh = 5000, directions = 8)
forclumpsSieve[is.na(forclumpsSieve)] <- 0

op <- par(mfrow=c(1, 2))
plot(rfForestMap, legend=FALSE)
plot(forclumpsSieve, legend=FALSE)
par(op)


# Make an NA-value raster based on the LC raster attributes
nonformask <- setValues(raster(rfForestMap), NA)
# Assign 0 to all cells corresponding to the non-forest class
nonformask[forclumpsSieve==0] <- 0

nonforclumpsSieve <- areaSieve(nonformask,thresh = 5000, directions = 8, keepzeros = TRUE)
nonforclumpsSieve[is.na(nonforclumpsSieve)] <- 1
plot(nonforclumpsSieve)

RFfinal <- nonforclumpsSieve
writeRaster(nonforclumpsSieve, filename = "C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map.tif", format= 'GTiff', overwrite=TRUE)



# This line was used to copy the extent of the SAR image in WGS84 to be later used in Google Earth Engine
library(rgdal)
SAR_wgs84 <- projectRaster(svmPredict, crs = crs('+init=epsg:4326'), method = "ngb")
SAR_wgs84


#### Validation of the FNF maps ####

# Loading validation points manually digitalized using ArcGIS Pro and Sentinel-2 imagery

crs(SAR_stack) <- crs('+init=epsg:31981')
ValidationPts <- read.csv("C:/Users/HP/Documents/Thesis/Data/Validation/validation_points_tile1291_table.csv", header=TRUE, sep = ",")
coordinates(ValidationPts)=~POINT_X+POINT_Y 
proj4string(ValidationPts)<-CRS("+init=epsg:31981")
ValidationPts@data <- ValidationPts@data[2]
ValidationPts@data[2:15] <- extract(SAR_stack, ValidationPts)
names(ValidationPts)[2:15] <-  c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh')
ValidationPts$class = as.factor(ValidationPts$class)

# Overall accuracy based on the comparison of the predicted class values and the class values manually assigned in ArcGiS
# Check every point with the basemap visually (if the classes are correct) !!!!!

RFfinal <- raster("C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map.tif")
crs(RFfinal) <- crs('+init=epsg:31981')
predicted = extract(RFfinal, ValidationPts)
predicted = as.factor(predicted)
confusionMatrix(predicted, ValidationPts$class) # Overall accuracy : 0.94

#  Compare the pixel-based map with object-based map (small study plot, more validation points in this area have to be added - 10 for non-forest and 40 for forest class)

ValidationPts <- read.csv("C:/Users/HP/Documents/Thesis/Data/Validation/validation_points_small_extent_table.csv", header=TRUE, sep = ",")
coordinates(ValidationPts)=~POINT_X+POINT_Y 
proj4string(ValidationPts)<-CRS("+init=epsg:31981")
ValidationPts@data <- ValidationPts@data[2]
ValidationPts@data[2:15] <- extract(SAR_stack, ValidationPts)
names(ValidationPts)[2:15] <-  c('avg_vv', 'md_vv', 'sd_vv', 'q10_vv', 'q90_vv', 'avg_f_vv', 'sd_f_vv', 'avg_vh', 'md_vh', 'sd_vh', 'q10_vh', 'q90_vh', 'avg_f_vh', 'sd_f_vh')
ValidationPts$class = as.factor(ValidationPts$class)

SegmFinal <- raster("C:/Users/HP/Documents/Thesis/Data/object_based_FNF_map_small.tif")
SegmFinal[SegmFinal == 2] <- 0

RFfinalCrop <- crop(RFfinal, SegmFinal)
writeRaster(RFfinalCrop, filename = "C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map_small.tif", format= 'GTiff', overwrite=TRUE)


predicted = extract(RFfinalCrop, ValidationPts)
predicted = as.factor(predicted)
confusionMatrix(predicted, ValidationPts$class)

predicted = extract(SegmFinal, ValidationPts)
predicted = as.factor(predicted)
confusionMatrix(predicted, ValidationPts$class)



#### Assessment using independent maps ####

# Load the CSV files with patch metrics derived from the predicted and independent maps using Fragstats software
# The Hansens map seems to be outdated (South America year 2000 fragments - maybe we should use another map?)

setwd('C:/Users/HP/Documents/Thesis/Data/Global FNF Maps')
# take the quanitile 95% or specify the limits of the histogram
FNFMapMetrics <- read.csv("FNFMapPatchMetrics.csv", header=TRUE, sep = ",")
TanDEMMetrics <- read.csv("TanDEMPatchMetrics.csv", header=TRUE, sep = ",")
MapbiomassMetrics <- read.csv("MapbiomassPatchMetrics.csv", header=TRUE, sep = ",")
HansensMetrics <- read.csv("HansensPatchMetrics.csv", header=TRUE, sep = ",")
JaxaMetrics <- read.csv("JaxaAlosPatchMetrics.csv", header=TRUE, sep = ",")
FNF30mMetrics <- read.csv("FNF30mMapPatchMetrics.csv", header=TRUE, sep = ",")

head(FNFMapMetrics)
head(TanDEMMetrics)
head(MapbiomassMetrics)
head(HansensMetrics)
head(JaxaMetrics)
head(FNF30mMetrics)

# Calculate the basic AREA statistics of all fragments
min(FNFMapMetrics$AREA)
max(FNFMapMetrics$AREA)
median(FNFMapMetrics$AREA)
mad(FNFMapMetrics$AREA)
quantile(FNFMapMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(TanDEMMetrics$AREA)
max(TanDEMMetrics$AREA)
median(TanDEMMetrics$AREA)
mad(TanDEMMetrics$AREA)
quantile(TanDEMMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(MapbiomassMetrics$AREA)
max(MapbiomassMetrics$AREA)
median(MapbiomassMetrics$AREA)
mad(MapbiomassMetrics$AREA)
quantile(MapbiomassMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(HansensMetrics$AREA)
max(HansensMetrics$AREA)
median(HansensMetrics$AREA)
mad(HansensMetrics$AREA)
quantile(HansensMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(JaxaMetrics$AREA)
max(JaxaMetrics$AREA)
median(JaxaMetrics$AREA)
mad(JaxaMetrics$AREA)
quantile(JaxaMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(FNF30mMetrics$AREA)
max(FNF30mMetrics$AREA)
median(FNF30mMetrics$AREA)
mad(FNF30mMetrics$AREA)
quantile(FNF30mMetrics$AREA, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)


# Calculate the basic SHAPE statistics of all fragments

min(FNFMapMetrics$FRAC)
max(FNFMapMetrics$FRAC)
median(FNFMapMetrics$FRAC)
mad(FNFMapMetrics$FRAC)
quantile(FNFMapMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(TanDEMMetrics$FRAC)
max(TanDEMMetrics$FRAC)
median(TanDEMMetrics$FRAC)
mad(TanDEMMetrics$FRAC)
quantile(TanDEMMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(MapbiomassMetrics$FRAC)
max(MapbiomassMetrics$FRAC)
median(MapbiomassMetrics$FRAC)
mad(MapbiomassMetrics$FRAC)
quantile(MapbiomassMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(HansensMetrics$FRAC)
max(HansensMetrics$FRAC)
median(HansensMetrics$FRAC)
mad(HansensMetrics$FRAC)
quantile(HansensMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(JaxaMetrics$FRAC)
max(JaxaMetrics$FRAC)
median(JaxaMetrics$FRAC)
mad(JaxaMetrics$FRAC)
quantile(JaxaMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)

min(FNF30mMetrics$FRAC)
max(FNF30mMetrics$FRAC)
median(FNF30mMetrics$FRAC)
mad(FNF30mMetrics$FRAC)
quantile(FNF30mMetrics$FRAC, probs = c(0.10, 0.25, 0.75, 0.99),na.rm=TRUE)



# The X and Y axes should be the same for all the histograms
# Histograms for the patch areas
hist(FNFMapMetrics[4:nrow(FNFMapMetrics), ]$AREA, breaks = function(x) {seq(from=0, to=5, by=0.5)}, main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(TanDEMMetrics[4:nrow(TanDEMMetrics), ]$AREA, breaks = function(x) {seq(from=0, to=20, by=0.5)}, main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(MapbiomassMetrics[4:nrow(MapbiomassMetrics), ]$AREA, breaks =  function(x) {seq(from=0, to=10, by=0.5)}, main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(HansensMetrics[4:nrow(HansensMetrics), ]$AREA, breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50), main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))



hist(FNFMapMetrics[2:nrow(FNFMapMetrics), ]$AREA, breaks = function(x) {seq(from=0, to=250, by=5)}, main='Distribution of patch areas', xlim=c(0, 250), xlab='Area [m^2]', ylim=c(0, 200))
hist(TanDEMMetrics[2:nrow(TanDEMMetrics), ]$AREA, breaks = {seq(from=0, to=250, by=5)}, main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(MapbiomassMetrics[2:nrow(MapbiomassMetrics), ]$AREA, breaks = {seq(from=0, to=250, by=5)}, main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(HansensMetrics[2:nrow(HansensMetrics), ]$AREA, breaks = {seq(from=0, to=250, by=5)}, main='Distribution of patch areas', 
     xlab='Area [m^2]', ylim=c(0, 200)) # remove the single forest clumps


# Histograms for the patch Fractal Dimension Index
hist(FNFMapMetrics$FRAC, main='Distribution of patch Fractal Dimension Index', breaks = c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35), xlab='FDI', ylim=c(0, 100))
hist(TanDEMMetrics$FRAC, main='Distribution of patch Fractal Dimension Index', breaks = c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35), xlab='FDI', ylim=c(0, 100))
hist(MapbiomassMetrics$FRAC, main='Distribution of patch Fractal Dimension Index', breaks = c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35), xlab='FDI', ylim=c(0, 100))
hist(HansensMetrics$FRAC, main='Distribution of patch Fractal Dimension Index', breaks = c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35), 
     xlab='FDI', ylim=c(0, 100)) # plot with the single pixels removed


# Aggregate my map into 30m resolution to get rid of the complex fragment shapes
FNF30m <- aggregate(RFfinal, fact = 3, fun = modal, na.rm = TRUE)
writeRaster(FNF30m, filename = "C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map_30m_res.tif", format= 'GTiff', overwrite = TRUE)

# Histograms for my map with pixels aggggregated to 30m resolution
hist(FNF30mMetrics$FRAC, main='Distribution of patch Fractal Dimension Index', breaks = c(1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35), xlab='FDI', ylim=c(0, 100))
hist(FNF30mMetrics[3:nrow(FNF30mMetrics), ]$AREA, breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50), main='Distribution of patch areas', xlab='Area [m^2]', ylim=c(0, 200))
hist(FNF30mMetrics[2:nrow(FNF30mMetrics), ]$AREA, breaks = function(x) {seq(from=0, to=250, by=5)}, main='Distribution of patch areas', xlim=c(0, 250), xlab='Area [m^2]', ylim=c(0, 200))



# Table min, max, sigma med, median, median absolute dev, q10, q99, q25, q75

# Visualise the fragments

library(rasterVis)

m1 <- raster('C:/Users/HP/Documents/Thesis/Data/Global FNF Maps/Mapbiomass_no_clumps.tif')
rf<-as.factor(m1)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="Mapbiomas", col.regions=c("beige", "dark green"))

m2 <- raster('C:/Users/HP/Documents/Thesis/Data/Global FNF Maps/TanDEM_no_clumps.tif')
rf<-as.factor(m2)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="TanDEM-X", col.regions=c("beige", "dark green"))

m3 <- raster('C:/Users/HP/Documents/Thesis/Data/Global FNF Maps/Hansen_no_clumps.tif')
rf<-as.factor(m3)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="GFC", col.regions=c("beige", "dark green"))

m4 <- raster('C:/Users/HP/Documents/Thesis/Data/Global FNF Maps/jaxa_alos.tif')
rf<-as.factor(m4)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="Jaxa ALOS", col.regions=c("beige", "dark green"))

m5 <- raster('C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map.tif')
rf<-as.factor(m5)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="FNF map", col.regions=c("beige", "dark green"))

m6 <- raster("C:/Users/HP/Documents/Thesis/Data/pixel_based_FNF_map_30m_res.tif")
rf<-as.factor(m6)
tar<-levels(rf)[[1]]
tar[["Landcover"]]<-c("Non-forest", "Forest")
levels(rf)<-tar
levelplot(rf,  main="FNF map 30m resolution", col.regions=c("beige", "dark green"))



