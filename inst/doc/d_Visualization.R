## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)

## ----generating data for plotting---------------------------------------------
library(voluModel) # Because of course
library(ggplot2) # For fancy plotting
library(viridisLite) # For high-contrast plotting palettes
library(dplyr) # To filter data
library(terra) # Now being transitioned in
library(sf) # Now being transitioned in

# Load data
oxygenSmooth <- rast(system.file("extdata/oxygenSmooth.tif", 
                                 package='voluModel'))

occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))

# Temperature
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- vect(paste0(td, "/temperature/woa18_decav_t00mn01_cropped.shp"))

# Creating a SpatRaster vector
template <- centerPointRasterTemplate(temperature)
tempTerVal <- rasterize(x = temperature, y = template, field = names(temperature))

# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(tempTerVal) <- envtNames
temperature <- tempTerVal

# Oxygen processing
names(oxygenSmooth) <- names(temperature)

# Clean points ----
occurrences <- occs %>% dplyr::select(decimalLongitude, decimalLatitude, depth) %>%
  distinct() %>% filter(dplyr::between(depth, 1, 2000))

# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(names(temperature))
occurrences$index <- unlist(lapply(occurrences$depth, 
                                   FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occurrences$index)
downsampledOccs <- data.frame()
for(i in indices){
  tempPoints <- occurrences[occurrences$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]], verbose = FALSE)
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  downsampledOccs <- rbind(downsampledOccs, tempPoints)
}
occsWdata <- downsampledOccs[,c("decimalLatitude", "decimalLongitude", "depth")] 

# Extract data ----
occsWdata$temperature <- xyzSample(occs = occsWdata, temperature)
occsWdata$AOU <- xyzSample(occs = occsWdata, oxygenSmooth)
occsWdata <- occsWdata[complete.cases(occsWdata),]

# Land shapefile
land <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")[1]

# Study region
studyRegion <- marineBackground(occsWdata, buff = 1000000)

# Get limits
tempLims <- quantile(occsWdata$temperature,c(0, 1))
aouLims <- quantile(occsWdata$AOU,c(0, 1))

# Reclassify environmental bricks to presence/absence
temperaturePresence <- classify(temperature, 
                                rcl = matrix(c(-Inf,tempLims[[1]],0,
                                               tempLims[[1]], tempLims[[2]], 1,
                                               tempLims[[2]], Inf, 0),
                                             ncol = 3, byrow = TRUE))
AOUpresence <- classify(oxygenSmooth, 
                        rcl = matrix(c(-Inf, aouLims[[1]],0,
                                       aouLims[[1]], aouLims[[2]], 1,
                                       aouLims[[2]], Inf, 0), 
                                     ncol = 3, byrow = TRUE))

# Put it all together
envelopeModel3D <- temperaturePresence * AOUpresence
envelopeModel3D <- mask(crop(envelopeModel3D, studyRegion), 
                        mask = studyRegion)
names(envelopeModel3D) <- names(temperature)
rm(AOUpresence, downsampledOccs, occurrences, temperaturePresence, 
   tempPoints, aouLims, envtNames, i, indices, layerNames, tempLims)

## ----pointMap, warning=FALSE, message=FALSE, eval=TRUE------------------------
pointMap(occs = occs, land = land, landCol = "black", spName = "Steindachneria argentea", 
         ptSize = 2, ptCol = "orange")

## ----pointCompMap, warning=FALSE, message=FALSE-------------------------------
pointCompMap(occs1 = occs, occs1Col = "red", occs1Name = "Raw", 
             occs2 = occsWdata, occs2Col = "orange", occs2Name = "Clean",
             spName = "Steindachneria argentea", agreeCol = "purple",
             land = land, landCol = "black", ptSize = 2, verbose = FALSE)

## ----oneRasterPlot------------------------------------------------------------
oneRasterPlot(rast = temperature[[1]],
              land = land, title = "Sea Surface Temperature, WOA 2018",
              landCol = "black", n = 11, option = "mako",
              varName = "Temperature")

## ----rasterComp---------------------------------------------------------------
rasterComp(rast1 = envelopeModel3D[[1]], rast1Name = "Surface",
           rast2 = envelopeModel3D[[10]], rast2Name = "45m", 
           land = land, landCol = "black", 
           title = "Suitability of Habitat for Luminous Hake\nAt Two Different Depths")

## ----plotLayers---------------------------------------------------------------
layerNames <- as.numeric(names(envelopeModel3D))
occsWdata$index <- unlist(lapply(occsWdata$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occsWdata$index)

layerPlot <- plotLayers(envelopeModel3D[[min(indices):max(indices)]], 
                        title = "Envelope Model of Luminous Hake,\n 20 to 700m",
                        land = land, landCol = "black")

## ----cleanup temporary directory----------------------------------------------
unlink(td, recursive = T)

