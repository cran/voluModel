## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 100)

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(voluModel) # Because of course
library(ggplot2) # For fancy plotting
library(rangeBuilder) # To compare marineBackground to getDynamicAlphaHull
library(dplyr) # To filter data
library(terra) # Now being transitioned in
library(sf) # Now being transitioned in

## ----occurrence dataset, message=FALSE, warning = FALSE-----------------------
# Get points
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))
occurrences <- occs %>% 
  dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  dplyr::distinct() %>% 
  dplyr::filter(depth %in% 1:2000)

## ----environmental data loading, eval=T, asis=T, message = FALSE, warning=FALSE----
# Temperature
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- vect(paste0(td, "/temperature/woa18_decav_t00mn01_cropped.shp"))

# Looking at the dataset
as.data.frame(temperature[1:5,1:10])

## ----temperature to compatible RasterBrick------------------------------------
# Creating a SpatRaster vector
template <- centerPointRasterTemplate(temperature)
tempTerVal <- rasterize(x = temperature, y = template, field = names(temperature))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(tempTerVal) <- envtNames
temperature <- tempTerVal

rm(tempTerVal)

## ----column interpretations, message=TRUE, warning=TRUE-----------------------
occsTest <- occurrences[19:24,]
xyzSample(occs = occsTest, envBrick = temperature)

colnames(occsTest) <- c("x", "y", "z")
xyzSample(occs = occsTest, envBrick = temperature)

rm(occsTest)

## ----downsample to voxel, eval=TRUE, warning=FALSE, message=FALSE-------------
# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(names(temperature))
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

## ----plot downsample, warning=FALSE-------------------------------------------
land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[1]

pointCompMap(occs1 = occs, occs2 = occurrences, 
             occs1Name = "original", occs2Name = "cleaned", 
             spName = "Steindachneria argentea", 
             land = land, verbose = FALSE)

## ----temperature extraction---------------------------------------------------
# Extraction
occurrences$temperature <- xyzSample(occs = occurrences, envBrick = temperature)

# Add "response" column for modeling
occurrences$response <- rep(1, times = nrow(occurrences))
occurrences <- occurrences[complete.cases(occurrences),]

head(occurrences)

## ----alpha hull demonstration, message=FALSE, warning=FALSE, eval=T-----------
trainingRegion <- marineBackground(occurrences, 
                                   fraction = 1, 
                                   partCount = 1, 
                                   buff = 1000000,
                                   clipToOcean = F)
plot(trainingRegion, border = F, col = "gray",
     main = "100% Points,\nMax 1 Polygon Permitted, 100 km Buffer",
     axes = T)
plot(land, col = "black", add = T)
points(occurrences[,c("decimalLongitude", "decimalLatitude")], 
       pch = 20, col = "red", cex = 1.5)

## ----clipToOcean demo, message=FALSE, warning=FALSE, eval=T-------------------
trainingRegion <- marineBackground(occurrences, 
                                   buff = 1000000, 
                                   clipToOcean = T)
plot(trainingRegion, border = F, col = "gray",
     main = "100 km Buffer,\nClipped to Occupied Polygon",
     axes = T)
plot(land, col = "black", add = T)
points(occurrences[,c("decimalLongitude", "decimalLatitude")], 
       pch = 20, col = "red", cex = 1.5)

## ----meridian wrap demo, warning=FALSE, message=FALSE, eval=T-----------------
# Fictional example occurrences
pacificOccs <- occurrences
pacificOccs$decimalLongitude <- pacificOccs$decimalLongitude - 100
for (i in 1:length(pacificOccs$decimalLongitude)){
  if (pacificOccs$decimalLongitude[[i]] < -180){
    pacificOccs$decimalLongitude[[i]] <- pacificOccs$decimalLongitude[[i]] + 360
  }
}

# marine Background
pacificTrainingRegion <- marineBackground(pacificOccs, 
                                          fraction = 0.95, partCount = 3,
                                          buff = 1000000,
                                          clipToOcean = T)
plot(pacificTrainingRegion, border = F, col = "gray",
     main = "marineBackground\nAntimeridian Wrap",
     axes = T)
plot(land, col = "black", add = T)
points(pacificOccs[,c("decimalLongitude", "decimalLatitude")], 
       pch = 20, col = "red", cex = 1.5)

## ----training points----------------------------------------------------------
# Background
backgroundVals <- mSampling3D(occs = occurrences, 
                              envBrick = temperature, 
                              mShp = trainingRegion, 
                              depthLimit = c(50, 1500))
backgroundVals$temperature <- xyzSample(occs = backgroundVals, temperature)

#Remove incomplete cases
backgroundVals <- backgroundVals[complete.cases(backgroundVals),]

## -----------------------------------------------------------------------------
# Add "response" column for modeling
backgroundVals$response <- rep(0, times = nrow(backgroundVals))

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
  scale_color_manual(values=c("#999999", "#E69F00"), 
                     labels = c("Pseudoabsence", "Presence")) +
  scale_fill_manual(values=c("#999999", "#E69F00"), 
                    labels = c("Pseudoabsence", "Presence")) +
  labs(title="Temperature Sampling Density,\nOccurrences vs. Pseudoabsences",
       x="Temperature (C)", y = "Density")+
  theme_classic()

## ----cleanup temporary directory----------------------------------------------
unlink(td, recursive = T)

