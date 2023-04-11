## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, error = FALSE, fig.retina = 1, dpi = 80)

## ----load packages, message=FALSE, warning=FALSE------------------------------
library(voluModel) # Since this is the package this vignette is about.
library(tibble) # For data organization
library(ggplot2) # For supplementary visualization
library(fields) # For raster interpolation
library(terra) # Now being transitioned in

## ----environmental data loading temperature, eval=T, message=FALSE, warning=FALSE, asis=T----
# Temperature
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- vect(paste0(td, "/temperature/woa18_decav_t00mn01_cropped.shp"))

# Looking at the dataset
head(temperature)

# Plotting the dataset
layout(matrix(c(1, 2), ncol=2, byrow=TRUE), widths=c(4, 1))
land <- rnaturalearth::ne_countries(scale = "small", 
                                    returnclass = "sf")[1]
temperatureForPlot <- temperature
crs(temperatureForPlot) <- crs(land) 
ext <- ext(temperatureForPlot)
plot(temperatureForPlot, main = "Distribution of voluModel Subset\nof WOA Temperature 2018",
     pch = 20, col = "red", xlim = ext[1:2], ylim = ext[3:4], cex = .6, mar = c(2,2,3,2))
plot(land, col = "black", add = T)

# What does the WOA depth structure look like?
depths <- names(temperatureForPlot)
depths <- as.numeric(gsub(depths[-1], pattern = "[d,M]", replacement = ""))
plot(0, xlim = c(0,1), ylim = c(0-max(depths), 0), axes=FALSE, type = "n", xlab = "", ylab = "Depth Intervals (m)")
axis(2, at = 0-depths, labels = depths)

## ----temperature processing, eval=FALSE---------------------------------------
#  # Creating a bottom raster
#  temperatureBottom <- bottomRaster(temperature)
#  
#  # Creating a SpatRaster vector
#  template <- centerPointRasterTemplate(temperature)
#  tempTerVal <- rasterize(x = temperature, y = template, field = names(temperature))
#  
#  # Get names of depths
#  envtNames <- gsub("[d,M]", "", names(temperature))
#  envtNames[[1]] <- "0"
#  names(tempTerVal) <- envtNames
#  temperature <- tempTerVal
#  rm(tempTerVal)
#  
#  # How do these files look?
#  par(mfrow=c(1,2))
#  p1 <- oneRasterPlot(temperature[[1]], land = land, landCol = "black",
#                title= "Surface Temperature (C)")
#  
#  p2 <- oneRasterPlot(temperatureBottom,land = land, landCol = "black",
#                title = "Bottom Temperature (C)")

## ----temperature plot, echo=FALSE, out.width = '100%', out.height= '100%'-----
knitr::include_graphics("TemperatureTopBottom.png")

## ----environmental data loading oxygen, eval=T, message=FALSE, warning=FALSE, asis=T----
td <- tempdir()
unzip(system.file("extdata/woa18_all_A00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/oxygen"), junkpaths = T)

oxygen <- vect(paste0(td, "/oxygen/woa18_all_A00mn01_cropped.shp"))

plot(oxygen, main = "Distribution of voluModel subset of WOA AOU 2018",
     pch = 20, col = "red", xlim = ext[1:2], ylim = ext[3:4], cex = .6)
plot(land, col = "black", add = T)

## ----interpolate oxygen, warning=FALSE, eval = F------------------------------
#  # Creating a SpatRaster vector
#  oxygen <- rasterize(x = oxygen, y = template,
#                     field = names(oxygen)) #Uses same raster template as temperature
#  
#  for (i in 1:nlyr(oxygen)){
#    oxygen[[i]] <- interpolateRaster(oxygen[[i]], lon.lat = T, fast = T, aRange = 5) #Thin plate spline interpolation
#    oxygen[[i]] <- crop(mask(x = oxygen[[i]],
#                             mask = temperature[[i]]),
#                        temperature[[i]])
#  }
#  
#  # Change names to match tempT
#  names(oxygen) <- envtNames

## ----smoothing oxygen, eval=FALSE---------------------------------------------
#  # Smoothing tempO and saving
#  oxygenSmooth <- oxygen
#  for (i in 1:nlyr(oxygenSmooth)){
#    oxygenSmooth[[i]] <- smoothRaster(oxygenSmooth[[i]], lon.lat = T) #Thin plate spline interpolation
#    oxygenSmooth[[i]] <- crop(mask(x = oxygenSmooth[[i]], mask = temperature[[i]]), temperature[[i]])
#  }
#  
#  # Change names to match tempT and save
#  names(oxygenSmooth) <- names(oxygen)
#  oxygenSmooth <- oxygenSmooth
#  
#  par(mfrow=c(1,2))
#  p3 <- oneRasterPlot(oxygen[[1]], land = land, landCol = "black",
#                title= "Surface Apparent Oxygen Utilization (µmol/kg),\ninterpolated")
#  p4 <- oneRasterPlot(oxygenSmooth[[1]], land = land, landCol = "black",
#       title = "Bottom Apparent Oxygen Utilization (µmol/kg),\ninterpolated and smoothed")

## ----AOU plot, echo=FALSE, out.width = '100%', out.height= '100%'-------------
knitr::include_graphics("AOUInterpAndSmooth.png")

## ----cleanup temporary directory----------------------------------------------
unlink(td, recursive = T)

