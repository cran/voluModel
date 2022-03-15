## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)
knitr::opts_knit$set(root.dir = system.file('extdata', 
                                            package='voluModel'))
load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))
load(system.file("extdata/oxygenBottom.RData",
                 package = 'voluModel'))

## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(voluModel) # Because of course
library(ggplot2) # For fancy plotting
library(rgdal) # For vector stuff. Will eventually be replaced with sf.
library(raster) # For raster stuff. Will eventually be replaced with terra.
library(viridisLite) # For high-contrast plotting palettes
library(latticeExtra) # Some fancy plotting

## ----loading data, message=FALSE, warning=FALSE, include=FALSE----------------
# Temperature
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- readOGR(dsn = paste0(td, "/temperature"), 
                       layer ="woa18_decav_t00mn01_cropped")
unlink(paste0(td, "/temperature"), recursive = T)

temperature@data[temperature@data == -999.999] <- NA

# Creating a bottom raster
temperatureBottom <- bottomRaster(temperature)

# Creating a RasterBrick
temperature <- rasterFromXYZ(cbind(temperature@coords,
                                   temperature@data))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(temperature) <- envtNames

# Oxygen processing
names(oxygenSmooth) <- names(temperature)

# Points
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))

#Clean points
occsClean <- occs[complete.cases(occs$depth),]
occsClean <- occsClean[occsClean$depth > 0.0,]
occsClean <- occsClean[occsClean$depth < 2000.0,]

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
occurrences <- occurrences[,c("decimalLatitude", "decimalLongitude", "depth")] 

## ----show points, warning=FALSE, echo=FALSE, message=FALSE--------------------
head(occurrences)

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
pointMap(occs = occurrences, ptCol = "orange", landCol = "black",
             spName = "Steindachneria argentia", ptSize = 3,
             land = land)
rm(indices, layerNames, tempPoints, i, downsampledOccs, occs)

## ----environmental data plotting, echo=FALSE, message=FALSE, warning=FALSE----
p1 <- oneRasterPlot(temperature[[1]], land = land, landCol = "black", 
              title= "Temperature (C)")
p2 <- oneRasterPlot(temperatureBottom,land = land, landCol = "black", 
              title = "Temperature (C)")

p3 <- oneRasterPlot(oxygenSmooth[[1]], land = land, landCol = "black", 
              title= "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")
p4 <- oneRasterPlot(oxygenBottom, land = land, landCol = "black",
     title = "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")

temp <- c("Surface" = p1, "Bottom" = p2)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                             horizontal = TRUE,
                             bg = "black",
                             fg = "white",
                             par.strip.text = list(col = "white", cex = 1.2, font = 2)))
temp <- c("Surface" = p3, "Bottom" = p4)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                             horizontal = TRUE,
                             bg = "black",
                             fg = "white",
                             par.strip.text = list(col = "white", cex = 1.2, font = 2)))

## ----data extraction, echo=FALSE, message=FALSE, warning=FALSE----------------
# Extract temperature
threeDimTemp <- xyzSample(occs = occurrences, temperature)
threeDimTemp <- cbind(rep("3D", length(threeDimTemp)), threeDimTemp)
surfTemp <- extract(x = temperature[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])
surfTemp <- cbind(rep("Surface", length(surfTemp)), surfTemp)
bottomTemp <- extract(x = temperatureBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])
bottomTemp <- cbind(rep("Bottom", length(bottomTemp)), bottomTemp)

# Extract AOU
threeDimOxy <- xyzSample(occs = occurrences, oxygenSmooth)
threeDimOxy <- cbind(rep("3D", length(threeDimOxy)), threeDimOxy)
surfOxy <- extract(x = oxygenSmooth[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])
surfOxy <- cbind(rep("Surface", length(surfOxy)), surfOxy)
bottomOxy <- extract(x = oxygenBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])
bottomOxy <- cbind(rep("Bottom", length(bottomOxy)), bottomOxy)

## ----temperature violin plots, echo=FALSE, message=FALSE, warning=FALSE-------
# Collation
allTemp <- rbind(threeDimTemp, surfTemp, bottomTemp)
colnames(allTemp) <- c("Group", "Temperature")
allTemp <- as.data.frame(allTemp)
allTemp$Temperature <- as.numeric(allTemp$Temperature)

allOxy <- rbind(threeDimOxy, surfOxy, bottomOxy)
colnames(allOxy) <- c("Group", "Oxygen")
allOxy <- as.data.frame(allOxy)
allOxy$Oxygen <- as.numeric(allOxy$Oxygen)

# Plotting
groups <- c("Surface", "Bottom", "3D")
tempPlot <- ggplot(allTemp, aes(x=Group, y=Temperature)) + 
  geom_boxplot(fill="#9C179EFF", notch = T) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(), 
        text = element_text(family = "Arial"), 
        axis.text = element_text(size = rel(1.1))) +
  scale_x_discrete(limits = groups) +
  ylab("Temperature (C)")

oxyPlot <- ggplot(allOxy, aes(x=Group, y=Oxygen)) + 
  geom_boxplot(fill="#ED7953FF", notch = T) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_blank(), 
        text = element_text(family = "Arial"), 
        axis.text = element_text(size = rel(1.1))) +
  scale_x_discrete(limits = groups) +
  ylab("Apparent Oxygen Utilization (µmol/kg)")

gridExtra::grid.arrange(tempPlot, oxyPlot, nrow = 1)

## ----study region, message=FALSE, warning=FALSE-------------------------------
studyRegion <- marineBackground(occurrences)
plot(temperature[[1]], 
     main = "Points and background sampling plotted on surface temperature",
     col = viridis(n = 11, option = "mako"))
plot(studyRegion, add = T, border = "orange", lwd = 2)
points(occurrences[,c("decimalLongitude","decimalLatitude")], 
       cex = 1, pch = 20, col = "red")

## ----sampling-----------------------------------------------------------------
# Surface Presences
oxyVals <- extract(x = oxygenSmooth[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])
tempVals <- extract(x = temperature[[1]], occurrences[,c("decimalLongitude", "decimalLatitude")])
vals <- cbind(occurrences, oxyVals, tempVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "AOU", "Temperature")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdataSurface <- vals

# Bottom Presences
oxyVals <- extract(x = oxygenBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])
tempVals <- extract(x = temperatureBottom, occurrences[,c("decimalLongitude", "decimalLatitude")])
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
temperaturePresence <- reclassify(temperature[[1]], 
                                  rcl = c(-Inf,tempLims[[1]],0,
                                          tempLims[[1]], tempLims[[2]], 1,
                                          tempLims[[2]], Inf, 0))
AOUpresence <- reclassify(oxygenSmooth[[1]], 
                          rcl = c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0))

# Put it all together, surface
envelopeModelSurface <- temperaturePresence * AOUpresence
envelopeModelSurface <- mask(crop(envelopeModelSurface, studyRegion), 
                             mask = studyRegion)

# Get limits, bottom
tempLims <- quantile(occsWdataBottom$Temperature,c(0, 1))
aouLims <- quantile(occsWdataBottom$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence, bottom
temperaturePresence <- reclassify(temperatureBottom, 
                                  rcl = c(-Inf,tempLims[[1]],0,
                                          tempLims[[1]], tempLims[[2]], 1,
                                          tempLims[[2]], Inf, 0))
AOUpresence <- reclassify(oxygenBottom, 
                          rcl = c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0))

# Put it all together, bottom
envelopeModelBottom <- temperaturePresence * AOUpresence
envelopeModelBottom <- mask(crop(envelopeModelBottom, studyRegion), mask = studyRegion)

## ----comparing 2D maps--------------------------------------------------------
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
rasterComp(rast1 = envelopeModelSurface, rast2 = envelopeModelBottom, 
           rast1Name = "Surface", rast2Name = "Bottom", land = land, landCol = "black", 
           title = "Comparison between surface and bottom envelope models")

## ----generate 3D envelope model-----------------------------------------------
# Get limits
tempLims <- quantile(occsWdata3D$Temperature,c(0, 1))
aouLims <- quantile(occsWdata3D$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence
temperaturePresence <- reclassify(temperature, 
                                  rcl = c(-Inf,tempLims[[1]],0,
                                          tempLims[[1]], tempLims[[2]], 1,
                                          tempLims[[2]], Inf, 0))
AOUpresence <- reclassify(oxygenSmooth, 
                          rcl = c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0))

# Put it all together
envelopeModel3D <- temperaturePresence * AOUpresence
envelopeModel3D <- mask(crop(envelopeModel3D, studyRegion), 
                        mask = studyRegion)
names(envelopeModel3D) <- names(temperature)

## ----plot 3D envelope model---------------------------------------------------
# Get indices of model-relevant layers
layerNames <- as.numeric(gsub("[X]", "", names(envelopeModel3D)))
occurrences$index <- unlist(lapply(occurrences$depth, 
                                   FUN = function(x) 
                                     which.min(abs(layerNames - x))))
indices <- sort(unique(occurrences$index))

plotLayers(envelopeModel3D[[min(indices):max(indices)]], 
          land = land)

