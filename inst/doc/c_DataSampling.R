## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)
load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))
load(system.file("extdata/oxygenBottom.RData",
                 package = 'voluModel'))

# Cleaning occurrences
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))

# Clean points
occsClean <- occs[complete.cases(occs$depth),]
occsClean <- occsClean[occsClean$depth > 0.0,]
occsClean <- occsClean[occsClean$depth < 2000.0,]
occurrences <- occsClean[,c("decimalLatitude", "decimalLongitude", "depth")] 
occurrences <- dplyr::distinct(occurrences)
occurrences <- occurrences[complete.cases(occurrences),]

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(voluModel) # Because of course
library(ggplot2) # For fancy plotting
library(rgdal) # For vector stuff. Will eventually be replaced with sf.
library(raster) # For raster stuff. Will eventually be replaced with terra.
library(rangeBuilder) # To compare marineBackground to getDynamicAlphaHull

## ----occurrence dataset, message=FALSE, warning = FALSE-----------------------
# Map occurrences
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
pointMap(occs = occurrences, ptCol = "orange", landCol = "black",
             spName = "Steindachneria argentia", ptSize = 3,
             land = land)

## ----alpha hull demonstration, message=FALSE, warning=FALSE, eval=FALSE-------
#  trainingRegion <- marineBackground(occurrences,
#                                  fraction = 1, partCount = 1, buff = 1000000,
#                                  clipToOcean = F)
#  plot(occurrences[,c("decimalLongitude", "decimalLatitude")],
#       main = "Mimimum of 100% Points in Training Region,\nMaximum of 1 Polygon Permitted, 100 km Buffer",
#       xlab = "Longitude", ylab = "Latitude",
#       pch = 20, col = "red",
#       xlim = c(-105, -40), ylim = c(0, 45))
#  plot(land, col = "black", add = T)
#  plot(trainingRegion, add = T, border = "orange", lwd = 2)

## ----plot alpha hull demonstration, echo=FALSE--------------------------------
knitr::include_graphics("alphaHullDemonstration-1.png")

## ----basic marineBackground demonstration, message=FALSE, warning=FALSE, eval=FALSE----
#  trainingRegion <- marineBackground(occurrences,
#                                  fraction = 1, partCount = 1,
#                                  clipToOcean = F)
#  plot(occurrences[,c("decimalLongitude", "decimalLatitude")],
#       main = "Mimimum of 100% Points in Training Region,\nMaximum of 1 Polygon Permitted, Automatic Buffer",
#       xlab = "Longitude", ylab = "Latitude",
#       pch = 20, col = "red",
#       xlim = c(-105, -40), ylim = c(0, 45))
#  plot(land, col = "black", add = T)
#  plot(trainingRegion, add = T, border = "orange", lwd = 2)

## ----plot basic marineBackground demonstration, echo=FALSE--------------------
knitr::include_graphics("basicMarineBackgroundDemonstration-1.png")

## ----clipToOcean demo, message=FALSE, warning=FALSE, eval=F-------------------
#  trainingRegion <- marineBackground(occurrences,
#                                  fraction = .95, partCount = 2, buff = 1000000,
#                                  clipToOcean = T)
#  plot(occurrences[,c("decimalLongitude", "decimalLatitude")],
#       main = "Minimum of 95% Points in Training Region,\nMaximum of 2 Polygons Permitted, 100 km Buffer",
#       xlab = "Longitude", ylab = "Latitude",
#       pch = 20, col = "red",
#       xlim = c(-105, -40), ylim = c(0, 45))
#  plot(land, col = "black", add = T)
#  plot(trainingRegion, add = T, border = "orange", lwd = 2)

## ----plot clipToOcean demo, echo=FALSE----------------------------------------
trainingRegion <- marineBackground(occurrences, 
                                fraction = .95, partCount = 2, buff = 1000000, 
                                clipToOcean = T)

knitr::include_graphics("clipToOceanDemo-1.png")

## ----meridian wrap demo, warning=FALSE, message=FALSE, eval=F-----------------
#  # Fictional example occurrences
#  pacificOccs <- occurrences
#  pacificOccs$decimalLongitude <- pacificOccs$decimalLongitude - 100
#  for (i in 1:length(pacificOccs$decimalLongitude)){
#    if (pacificOccs$decimalLongitude[[i]] < -180){
#      pacificOccs$decimalLongitude[[i]] <- pacificOccs$decimalLongitude[[i]] + 360
#    }
#  }
#  
#  # marine Background
#  pacificTrainingRegion <- marineBackground(pacificOccs,
#                                            fraction = 1, partCount = 2,
#                                            clipToOcean = T)
#  plot(pacificOccs[,c("decimalLongitude", "decimalLatitude")],
#       main = "marineBackground Antimeridian Wrap",
#       xlab = "Longitude", ylab = "Latitude",
#       ylim = c(-40,60),
#       pch = 20, col = "red")
#  plot(land, col = "black", add = T)
#  plot(pacificTrainingRegion, add = T, border = "orange", lwd = 2)

## ----plot meridian wrap demo, echo=FALSE--------------------------------------
knitr::include_graphics("meridianWrapDemo-1.png")

## ----environmental data loading, eval=T, asis=T, message = FALSE, warning=FALSE----
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- readOGR(dsn = paste0(td, "/temperature"), 
                       layer ="woa18_decav_t00mn01_cropped")
unlink(paste0(td, "/temperature"), recursive = T)

temperature@data[temperature@data == -999.999] <- NA

# How are the data structured?
head(temperature@data)

## ----temperature to compatible RasterBrick------------------------------------
# Creating a RasterBrick
temperature <- rasterFromXYZ(cbind(temperature@coords,
                                   temperature@data))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(temperature) <- envtNames

# Here's a sampling of depth plots from the 102 depth layers available
plot(temperature[[c(1, 10, 50, 100)]])

## ----column interpretations, message=TRUE, warning=TRUE-----------------------
occsTest <- occurrences[1:5,]
xyzSample(occs = occsTest, envBrick = temperature)
colnames(occsTest) <- c("x", "y", "z")
xyzSample(occs = occsTest, envBrick = temperature)

# Of course, the names have to be somewhat interpretable
colnames(occsTest) <- c("ham", "eggs", "cheese")
xyzSample(occs = occsTest, envBrick = temperature)
rm(occsTest)

## ----downsample to voxel, eval=TRUE, warning=FALSE, message=FALSE-------------
occurrences <- occsClean[,c("decimalLatitude", "decimalLongitude", "depth")] 

# Preliminary cleaning
occurrences <- dplyr::distinct(occurrences)
occurrences <- occurrences[complete.cases(occurrences),]

# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
occurrences$index <- unlist(lapply(occurrences$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occurrences$index)
downsampledOccs <- data.frame()
for(i in indices){
  tempPoints <- occurrences[occurrences$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  downsampledOccs <- rbind(downsampledOccs, tempPoints)
}

occurrences <- downsampledOccs
head(occurrences)

print(paste0("Original number of points: ", nrow(occs), "; number of downsampled occs: ", nrow(occurrences)))

## ----plot downsample----------------------------------------------------------
pointCompMap(occs1 = occs, occs2 = occurrences, 
             occs1Name = "Original", occs2Name = "Cleaned", 
             spName = "Steindachneria argentia", 
             land = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1])

## ----temperature extraction---------------------------------------------------
# Extraction
occurrences$temperature <- xyzSample(occs = occurrences, envBrick = temperature)

# Add "response" column for modeling
occurrences$response <- rep(1, times = nrow(occurrences))
occurrences <- occurrences[complete.cases(occurrences),]

head(occurrences)

## ----training points----------------------------------------------------------
# Background
backgroundVals <- mSampling3D(occs = occurrences, 
                              envBrick = temperature, 
                              mShp = trainingRegion, 
                              depthLimit = c(50, 1500))
backgroundVals$temperature <- xyzSample(occs = backgroundVals, temperature)

## -----------------------------------------------------------------------------
# Add "response" column for modeling
backgroundVals$response <- rep(0, times = nrow(backgroundVals))
backgroundVals <- backgroundVals[complete.cases(backgroundVals),]

head(backgroundVals)

## ----glm data prep part 2-----------------------------------------------------
# Sample background points weighted by distance from mean of occurrence temperatures
meanTemp <- mean(occurrences[,c("temperature")])
backgroundVals$distance <- abs(backgroundVals[,"temperature"] - meanTemp)
backgroundVals$sampleWeight <- (backgroundVals$distance - 
                                  min(backgroundVals$distance))/(max(backgroundVals$distance) -
                                                                  min(backgroundVals$distance))
sampleForAbsence <- sample(x = rownames(backgroundVals), 
                           size = nrow(occurrences) * 100, 
                           prob = backgroundVals$backgroundVals)
backgroundVals <- backgroundVals[match(sampleForAbsence, 
                                       rownames(backgroundVals)),]

backgroundVals$response <- as.factor(backgroundVals$response)

# Unite datasets and see how things look
dataForModeling <- rbind(occurrences, backgroundVals[,colnames(occurrences)])
ggplot(dataForModeling, aes(x = temperature, fill = response, color = response)) +
  geom_density(alpha = .6) +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                     labels = c("Pseudoabsence", "Presence")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                    labels = c("Pseudoabsence", "Presence")) +
  labs(title="Temperature Sampling Density,\nOccurrences vs. Pseudoabsences",
       x="Temperature (C)", y = "Density")+
  theme_classic()

