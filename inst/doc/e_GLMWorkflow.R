## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)
knitr::opts_knit$set(root.dir = system.file('extdata', 
                                            package='voluModel'))
load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))

## ----packages, warning=FALSE--------------------------------------------------
library(voluModel) # Because of course
library(dplyr) # For occurrence data filtration
library(ggplot2) # For fancy plotting
library(rgdal) # For vector stuff. Will eventually be replaced with sf.
library(raster) # For raster stuff. Will eventually be replaced with terra.
library(viridisLite) # For high-contrast plotting palettes

## ----occurrence data, eval=T--------------------------------------------------
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))
summary(occs$depth)
boxplot.stats(occs$depth)

## ----clean occurrence data, eval=T--------------------------------------------
occsClean <- occs %>% dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  distinct() %>% filter(depth %in% 1:2000)
head(occsClean)

ggplot(occsClean, aes(x = 1, y=depth)) +
  geom_violin(fill="yellow", ) +
  theme_classic(base_size = 15) +
  theme(axis.title = element_blank(),
        text = element_text(family = "Arial"),
        axis.text = element_text(size = rel(1.1))) +
  ggtitle("Depth") +
  geom_boxplot(width=0.1)

## ----environmental data loading, eval=T, asis=T, message=FALSE----------------
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- readOGR(dsn = paste0(td, "/temperature"), 
                       layer ="woa18_decav_t00mn01_cropped")
unlink(paste0(td, "/temperature"), recursive = T)

temperature@data[temperature@data == -999.999] <- NA

# Creating a RasterBrick
temperature <- rasterFromXYZ(cbind(temperature@coords,
                                   temperature@data))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(temperature) <- envtNames

## ----environmental data loading oxygen, eval=T, asis=T, warning=FALSE---------
load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))

# Change names to match temperature
names(oxygenSmooth) <- names(temperature)

## ----downsample to voxel, eval=TRUE, warning=FALSE, message=FALSE-------------
# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
occsClean$index <- unlist(lapply(occsClean$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occsClean$index)

# Downsamples occsClean in each depth layer
downsampledOccs <- data.frame()
for(i in indices){
  tempPoints <- occsClean[occsClean$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  downsampledOccs <- rbind(downsampledOccs, tempPoints)
}

occsClean <- downsampledOccs

print(paste0("Original number of points: ", nrow(occs), "; number of downsampled occs: ", nrow(occsClean)))
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
pointCompMap(occs1 = occs, occs2 = occsClean, 
             occs1Name = "Original", occs2Name = "Cleaned", 
             spName = "Steindachneria argentea", 
             land = land)

## ----environmental background sampling, warning=FALSE-------------------------
backgroundSamplingRegions <- marineBackground(occsClean, buff = 1000000)
plot(temperature[[1]], 
     main = "Points and background sampling plotted on surface temperature",
     col = viridis(n= 11, option = "mako"))
plot(backgroundSamplingRegions, add = T, border = "orange", lwd = 2)
points(occsClean[,c("decimalLongitude","decimalLatitude")], 
       cex = 1, pch = 20, col = "red")

## ----presence sampling--------------------------------------------------------
# Presences
oxyVals <- xyzSample(occs = occsClean, oxygenSmooth)
tempVals <- xyzSample(occs = occsClean, temperature)
vals <- cbind(occsClean, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdata <- vals

# Add response as a column
occsWdata$response <- rep(1, times = nrow(occsWdata))

## ----background sampling------------------------------------------------------
# Background
backgroundVals <- mSampling3D(occs = occsClean, 
                              envBrick = temperature, 
                              mShp = backgroundSamplingRegions, 
                              depthLimit = c(5, 800))
oxyVals <- xyzSample(occs = backgroundVals, oxygenSmooth)
tempVals <- xyzSample(occs = backgroundVals, temperature)
vals <- cbind(backgroundVals, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
backgroundWdata <- vals

# Add response as a column
backgroundWdata$response <- rep(0, times = nrow(backgroundWdata))

## ----glm data prep part 2-----------------------------------------------------
# Sample background points weighted by distance from centroid of occurrence environments
suitableCentroid <- apply(occsWdata[,c("Temperature", "AOU")], 
                          MARGIN = 2, FUN = mean)
backgroundWdata$distance <- apply(backgroundWdata[,c("Temperature", "AOU")], MARGIN = 1, 
                                  FUN = function(x) dist(rbind(suitableCentroid, x)))
backgroundWdata$sampleWeight <- (backgroundWdata$distance - 
                                 min(backgroundWdata$distance))/(max(backgroundWdata$distance)-
                                                                 min(backgroundWdata$distance))
sampleForAbsence <- sample(x = rownames(backgroundWdata), 
                           size = nrow(occsWdata) * 100, 
                           prob = backgroundWdata$sampleWeight)
backgroundWdata <- backgroundWdata[match(sampleForAbsence, 
                                         rownames(backgroundWdata)),]

## ----uniting datasets---------------------------------------------------------
# Unite datasets
datForMod <- rbind(occsWdata, backgroundWdata[,colnames(occsWdata)])
rm(suitableCentroid, sampleForAbsence, backgroundWdata, occsWdata)

## ----generate glm niche model-------------------------------------------------
glmModel <- glm(formula = response ~ Temperature *  AOU, 
                  family = binomial(link = "logit"),  data = datForMod)
summary(glmModel)

## ----project glm niche model--------------------------------------------------
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
index <- seq(from = match(min(datForMod$depth), layerNames), 
             to = match(max(datForMod$depth), layerNames), by = 1)
depthPred <- NULL
for(j in index){
  depthPreds <- stack(temperature[[j]], oxygenSmooth[[j]])
  crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(depthPreds) <- c("Temperature", "AOU")
  depthPred[[j]] <- mask(predict(depthPreds, glmModel), backgroundSamplingRegions)
  depthPred[[j]] <- crop(depthPred[[j]], backgroundSamplingRegions)
  names(depthPred[[j]]) <- layerNames[[j]]
}
glmPred <- stack(depthPred[!unlist(lapply(depthPred, FUN = function(x) is.null(x)))])

## ----glm niche threshold------------------------------------------------------
glmThreshold <- quantile(xyzSample(datForMod[datForMod$response == 1,
                                             c("depth", "decimalLongitude", "decimalLatitude")], 
                                   brick(glmPred)), .1, na.rm = T)[[1]] # MS90
glmThresholded <- glmPred > glmThreshold
glmThresholded <- reclassify(glmThresholded, c(NA, NA, 0), include.lowest = T)

## ----thresholded glm niche model plotted--------------------------------------
plotLayers(glmThresholded, 
          land = land, landCol = "black",
          title = "Areas of suitable habitat for \n Luminous Hake, 5m to 800m")

## ----calculate MESS, warning=FALSE--------------------------------------------
# Prepare environmental data
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
datForMod$index <- unlist(lapply(datForMod$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(datForMod$index)

projList <- list("AOU" = oxygenSmooth[[min(indices):max(indices)]], 
                 "Temperature" = temperature[[min(indices):max(indices)]])

# Calculate MESS
messBrick <- MESS3D(calibration = datForMod, projection = projList)

## ----reclassify and plot extrapolation----------------------------------------
# Reclassify MESS 
extrapolation <- 0 >= messBrick

# Plot Extrapolation
plotLayers(extrapolation, land = land, landCol = "black", 
           title = "Areas of extrapolation according to MESS,\n 5m to 800m")

## ----reclassify MESS and plot GLM with no extrapolation, warning=FALSE--------
# Reclassify MESS 
noExtrapolation <- messBrick > 0

glmThreshNoExtrapolation <- NULL
for(i in 1:nlayers(noExtrapolation)){
  croppedLayer <- glmThresholded[[i]] * noExtrapolation[[i]]
  glmThreshNoExtrapolation <- stack(c(glmThreshNoExtrapolation, croppedLayer))
}

# Plot MESS
plotLayers(glmThreshNoExtrapolation, land = land, 
           landCol = "black", 
           title = "Areas of suitable habitat with no model extrapolation,\n 5m to 800m")

