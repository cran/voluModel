## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(voluModel) # Because of course
library(dplyr) # To filter data
library(ggplot2) # For fancy plotting
library(terra) # Now being transitioned in

## ----show points, warning=FALSE, echo=TRUE, message=FALSE, eval=TRUE----------
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))

# Filter points
occurrences <- occs %>% 
  dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  dplyr::distinct() %>% 
  filter(depth %in% 1:2000)

head(occurrences)
land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[1]
pointMap(occs = occurrences, ptCol = "orange", landCol = "black",
             spName = "Steindachneria argentea", ptSize = 3,
             land = land)

## ----loading temperature data, message=FALSE, warning=FALSE, include=TRUE-----
# Temperature
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- vect(paste0(td, "/temperature/woa18_decav_t00mn01_cropped.shp"))

# Creating a bottom raster
temperatureBottom <- bottomRaster(temperature)

# Creating a SpatRaster vector
template <- centerPointRasterTemplate(temperature)
tempTerVal <- terra::rasterize(x = temperature, y = template, field = names(temperature))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(tempTerVal) <- envtNames
temperature <- tempTerVal
rm(tempTerVal)

## ----plotting temperature, eval = FALSE---------------------------------------
#  # How do these files look?
#  par(mfrow=c(1,2))
#  p1 <- oneRasterPlot(temperature[[1]], land = land, landCol = "black",
#                title= "Surface Temperature (C)")
#  
#  p2 <- oneRasterPlot(temperatureBottom,land = land, landCol = "black",
#                title = "Bottom Temperature (C)")

## ----temperature plot, echo=FALSE, out.width = '100%', out.height= '100%'-----
knitr::include_graphics("TemperatureTopBottom.png")

## ----loading oxygen data, message=FALSE, warning=FALSE, include=TRUE----------
# Oxygen processing, pre-baked to save time
oxygenSmooth <- rast(system.file("extdata/oxygenSmooth.tif", 
                                 package='voluModel'))
oxygenBottom <- rast(system.file("extdata/oxygenBottom.tif", 
                                 package="voluModel"))
names(oxygenSmooth) <- names(temperature)

## ----oxygen plotting, echo=TRUE, eval = FALSE, message=FALSE, warning=FALSE----
#  par(mfrow=c(1,2))
#  p3 <- oneRasterPlot(oxygenSmooth[[1]], land = land, landCol = "black",
#                title= "Surface Apparent Oxygen Utilization (µmol/kg),\ninterpolated and smoothed")
#  p4 <- oneRasterPlot(oxygenBottom, land = land, landCol = "black",
#       title = "Bottom Apparent Oxygen Utilization (µmol/kg),\ninterpolated and smoothed")

## ----oxygen plot, echo=FALSE, out.width = '100%', out.height= '100%'----------
knitr::include_graphics("OxygenTopBottom.png")

## ----occurrence and depth matchup---------------------------------------------
# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(names(temperature))
occurrences$index <- unlist(lapply(occurrences$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occurrences$index)

# Downsamples occurrences in each depth layer
downsampledOccs <- data.frame()
for(i in indices){
  tempPoints <- occurrences[occurrences$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]], verbose = FALSE)
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  downsampledOccs <- rbind(downsampledOccs, tempPoints)
}
occurrences <- downsampledOccs
occurrences <- occurrences[,c("decimalLatitude", "decimalLongitude", "depth")]

rm(indices, layerNames, tempPoints, i, downsampledOccs, occs)

## ----data extraction, echo=FALSE, message=FALSE, warning=FALSE----------------
# Extract temperature
threeDimTemp <- xyzSample(occs = occurrences, temperature)
threeDimTemp <- cbind(rep("X, Y,\nZ", length(threeDimTemp)), threeDimTemp)
surfTemp <- extract(x = temperature[[1]], 
                    occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
surfTemp <- cbind(rep("X, Y,\nSurface", length(surfTemp)), surfTemp)
bottomTemp <- extract(x = temperatureBottom, 
                      occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
bottomTemp <- cbind(rep("X, Y,\nBottom", length(bottomTemp)), bottomTemp)

# Extract AOU
threeDimOxy <- xyzSample(occs = occurrences, oxygenSmooth)
threeDimOxy <- cbind(rep("X, Y,\nZ", length(threeDimOxy)), threeDimOxy)
surfOxy <- extract(x = oxygenSmooth[[1]], 
                   occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
surfOxy <- cbind(rep("X, Y,\nSurface", length(surfOxy)), surfOxy)
bottomOxy <- extract(x = oxygenBottom, 
                     occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
bottomOxy <- cbind(rep("X, Y,\nBottom", length(bottomOxy)), bottomOxy)

## ----temperature violin plots, echo=FALSE, message=FALSE, warning=FALSE-------
# Collation
allTemp <- rbind(threeDimTemp, surfTemp, bottomTemp)
colnames(allTemp) <- c("Group", "Temperature")
allTemp <- as.data.frame(allTemp)
allTemp <- allTemp[complete.cases(allTemp),]
allTemp$Temperature <- as.numeric(allTemp$Temperature)

allOxy <- rbind(threeDimOxy, surfOxy, bottomOxy)
colnames(allOxy) <- c("Group", "Oxygen")
allOxy <- as.data.frame(allOxy)
allOxy <- allOxy[complete.cases(allOxy),]
allOxy$Oxygen <- as.numeric(allOxy$Oxygen)

# Plotting
groups <- c("X, Y,\nSurface", "X, Y,\nBottom", "X, Y,\nZ")
tempPlot <- ggplot(allTemp, aes(x=Group, y=Temperature)) + 
  geom_boxplot(fill="#b2182b", notch = TRUE) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(), 
        text = element_text(family = "Arial"), 
        axis.text = element_text(size = rel(1.1))) +
  scale_x_discrete(limits = groups) +
  ylab("Temperature (C)")

oxyPlot <- ggplot(allOxy, aes(x=Group, y=Oxygen)) + 
  geom_boxplot(fill="#2166ac", notch = TRUE) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(), 
        text = element_text(family = "Arial"), 
        axis.text = element_text(size = rel(1.1))) +
  scale_x_discrete(limits = groups) +
  ylab("Apparent Oxygen Utilization (µmol/kg)")

gridExtra::grid.arrange(tempPlot, oxyPlot, nrow = 1)

## ----study region, message=FALSE, warning=FALSE, eval=FALSE-------------------
#  studyRegion <- marineBackground(occurrences, buff = 1000000, clipToOcean = TRUE)
#  landVect <- vect(land)
#  landVect <- terra::project(landVect, y = studyRegion)
#  plot(studyRegion, border = F, col = "gray",
#       main = "Points and Background Sampling",
#       axes = T)
#  plot(landVect, col = "black", add = T)
#  points(occurrences[,c("decimalLongitude", "decimalLatitude")],
#         pch = 20, col = "red", cex = 1.5)

## ----study region hidden, message=FALSE, warning=FALSE, eval=TRUE, echo=FALSE----
studyRegion <- vect(system.file("extdata/backgroundSamplingRegions.shp",
                                package='voluModel'))

## ----plot study region, echo=FALSE, out.width = '100%', out.height= '100%'----
knitr::include_graphics("PointsAndTrainingRegion.png")

## ----sampling-----------------------------------------------------------------
# Surface Presences
oxyVals <- extract(x = oxygenSmooth[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
tempVals <- extract(x = temperature[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
vals <- cbind(occurrences, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdataSurface <- vals

# Bottom Presences
oxyVals <- extract(x = oxygenBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
tempVals <- extract(x = temperatureBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])[,2]
vals <- cbind(occurrences, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdataBottom <- vals

# 3D Presences
oxyVals <- xyzSample(occs = occurrences, oxygenSmooth)
tempVals <- xyzSample(occs = occurrences, temperature)
vals <- cbind(occurrences, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdata3D <- vals

rm(oxyVals, tempVals, vals)

## ----generate 2D envelope niche models----------------------------------------
# Get limits, surface
tempLims <- quantile(occsWdataSurface$Temperature,c(0, 1))
aouLims <- quantile(occsWdataSurface$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence, surface
temperaturePresence <- classify(temperature[[1]], 
                                  rcl = matrix(c(-Inf,tempLims[[1]],0, 
                                                 tempLims[[1]], tempLims[[2]], 1,
                                                 tempLims[[2]], Inf, 0), 
                                               ncol = 3, byrow = TRUE))
AOUpresence <- classify(oxygenSmooth[[1]], 
                          rcl = matrix(c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0), ncol = 3, byrow = TRUE))

# Put it all together, surface
envelopeModelSurface <- temperaturePresence * AOUpresence
envelopeModelSurface <- mask(crop(envelopeModelSurface, studyRegion), 
                             mask = studyRegion)

# Get limits, bottom
tempLims <- quantile(occsWdataBottom$Temperature,c(0, 1))
aouLims <- quantile(occsWdataBottom$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence, bottom
temperaturePresence <- classify(temperatureBottom, 
                                  rcl = matrix(c(-Inf,tempLims[[1]],0,
                                          tempLims[[1]], tempLims[[2]], 1,
                                          tempLims[[2]], Inf, 0), byrow = TRUE, ncol=3))
AOUpresence <- classify(oxygenBottom, 
                          rcl = matrix(c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0), byrow = TRUE, ncol = 3))

# Put it all together, bottom
envelopeModelBottom <- temperaturePresence * AOUpresence
envelopeModelBottom <- mask(crop(envelopeModelBottom, studyRegion), mask = studyRegion)

## ----comparing 2D maps--------------------------------------------------------
rasterComp(rast1 = envelopeModelSurface, rast2 = envelopeModelBottom, 
           rast1Name = "Surface", rast2Name = "Bottom", 
           land = land, landCol = "black", 
           title = "Comparison between surface and bottom envelope models")

## ----generate 3D envelope model-----------------------------------------------
# Get limits
tempLims <- quantile(occsWdata3D$Temperature,c(0, 1))
aouLims <- quantile(occsWdata3D$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence
temperaturePresence <- classify(temperature, 
                                  rcl = matrix(c(-Inf,tempLims[[1]],0,
                                                 tempLims[[1]], tempLims[[2]], 1,
                                                 tempLims[[2]], Inf, 0), 
                                               byrow = TRUE, ncol = 3))
AOUpresence <- classify(oxygenSmooth, 
                          rcl = matrix(c(-Inf, aouLims[[1]],0,
                                         aouLims[[1]], aouLims[[2]], 1,
                                         aouLims[[2]], Inf, 0), 
                                       byrow = TRUE, ncol = 3))

# Put it all together
envelopeModel3D <- temperaturePresence * AOUpresence
envelopeModel3D <- mask(crop(envelopeModel3D, studyRegion), 
                        mask = studyRegion)
names(envelopeModel3D) <- names(temperature)

## ----plot 3D envelope model---------------------------------------------------
# Get indices of model-relevant layers
layerNames <- as.numeric(names(envelopeModel3D))
occurrences$index <- unlist(lapply(occurrences$depth, 
                                   FUN = function(x) 
                                     which.min(abs(layerNames - x))))
indices <- sort(unique(occurrences$index))

plotLayers(envelopeModel3D[[min(indices):max(indices)]],
           title = "Envelope Model of Luminous Hake,\n 20 to 700m",
           land = land)

## ----cleanup temporary directory----------------------------------------------
unlink(td, recursive = T)

