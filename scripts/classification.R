library(tidyverse)
library(sf)
library(terra)
library(randomForest)
library(caret)

# Add path to the folder where PlanetScope imagery is stored
rast_dir <- ""

# Preprocess teaching data from the GHG plots

e <- read.csv("data/veg_class.csv", sep=";")
e18 <- read.csv("data/veg_class_extra18.csv")
e17 <- read.csv("data/veg_class_extra17.csv")

# Coordinates for the plots
c17_extra <- read.csv("data/veg_class_extra17_coords.csv")
c18 <- read.csv("data/coords_ghg_doc_co2_epsg3067.csv") %>% 
  mutate(plot = as.character(plot))
# Select the extra coordinates that we are using for 2018
c18_extra <- c18 %>% filter(c18$plot %in% e18$plot)

# Modify plot names
e18$plot <- paste(e18$plot, "extra", sep="_")
e17$plot <- paste(e17$plot, "extra", sep="_")

c18_extra$plot <- paste(c18_extra$plot, "extra", sep="_")
c17_extra$plot <- paste(c17_extra$plot, "extra", sep="_")

# Merge all three vegetation classification files
d <- rbind(e, e18, e17)

# Merge coordinate files so that the data frame only has plot, x, y, (in epsg 3067) and Comment columns for veg classification plots.
coords <- data.frame(matrix(ncol=4))
colnames(coords) <- c("plot", "x", "y", "Comments")

coords <- coords %>% add_row(plot=c18$plot, x=c18$Easting_GHG, y=c18$Northing_GHG, Comments=paste(c18$Comment_GHG))
coords <- coords %>% add_row(plot=c18_extra$plot, x=c18_extra$Easting_CO2, y=c18_extra$Northing_CO2, Comments=paste(c18_extra$Comment_CO2))
coords <- coords %>% add_row(plot=c17_extra$plot, x=c17_extra$x, y=c17_extra$y, Comments=NA)

# Remove columns that don't have coordinate data
coords <- subset(coords, !is.na(x))

# Merge coords with veg
coords_veg <- merge(coords, d, by="plot") 
final <- subset(coords_veg, select=c(plot, cavm, cavm.notes, x, y, Comments))

final %>% mutate(cl = as.numeric(as.factor(cavm))) -> d

# GHG dataset to spatial points
p <- st_as_sf(d, coords = c("x", "y"), crs = 3067)


# Veg classifications from a previous veg survey
c16 <- read.csv2("data/2x2habitaatit.csv", stringsAsFactors = F, dec = ".")

# Transform to long format, create plot-level coordinates and filter relevant classes
# Recode the classes to match the classes in the GHG dataset
c16 %>% 
  pivot_longer(cols = C:W, names_to = "Plot", values_to = "class") %>% 
  mutate(x = ifelse(Plot == "E", x + 5.5, x)) %>% 
  mutate(x = ifelse(Plot == "W", x - 5.5, x)) %>% 
  mutate(y = ifelse(Plot == "N", y + 5.5, y)) %>% 
  mutate(y = ifelse(Plot == "S", y - 5.5, y)) %>% 
  mutate(site = ID) %>% 
  mutate(ID = paste0(ID, Plot)) %>% 
  arrange(ID) %>% 
  filter(class %in% c(12,22,23,24,25,27,29,31,32,80,81,85,88,89)) %>% 
  mutate(cavm = as.character(class)) %>% 
  mutate(cavm = recode(cavm, "12" = "es",
                       "22" = "es",
                       "23" = "ds",
                       "24" = "ds",
                       "25" = "es",
                       "27" = "es",
                       "29" = "g",
                       "31" = "g",
                       "32" = "g",
                       "80" = "g",
                       "81" = "b",
                       "85" = "g",
                       "88" = "w",
                       "89" = "w")) -> p2

# Add the final classifications conducted by interpreting aerial imagery 
p3 <- st_read("data/added_vegpoints.gpkg") %>% 
  filter(!is.na(class)) %>% 
  st_transform(crs = 32634)

p3 <- p3 %>% 
  mutate(class = recode(class, "1" = "b",
                       "2" = "es",
                       "3" = "ds",
                       "4" = "g",
                       "5" = "w")) %>% 
  rename(cavm = class)


#####################################################
# Preprocess satellite imagery

# The extent of the study area
aoi <- st_read("data/larger_landscapewindow_extent.shp")

# The imagery from September consists of two scenes
# clip and combine together ->
r31 <- rast(paste0(rast_dir, "20180909_094848_1027_3B_AnalyticMS_SR.tif"))
r32 <- rast(paste0(rast_dir, "20180909_094849_1027_3B_AnalyticMS_SR.tif"))

aoi <- st_transform(aoi, crs = crs(r31, proj = T))
p <- st_transform(p, crs = crs(r31, proj = T))

r31 <- crop(r31, aoi)
r32 <- crop(r32, aoi)

r3 <- mosaic(r31, r32, fun = "mean")

# Clip and combine the early July imagery
r11 <- rast(paste0(rast_dir, "20180702_094658_102c_3B_AnalyticMS_SR.tif"))
r12 <- rast(paste0(rast_dir, "20180702_094657_102c_3B_AnalyticMS_SR.tif"))

r11 <- crop(r11, aoi)
r12 <- crop(r12, aoi)

r1 <- mosaic(r11, r12, fun = "mean")

# Clip and combine the late July imagery
r21 <- rast(paste0(rast_dir, "20180729_094920_1006_3B_AnalyticMS_SR.tif"))
r22 <- rast(paste0(rast_dir, "20180729_094921_1006_3B_AnalyticMS_SR.tif"))

r21 <- crop(r21, aoi)
r22 <- crop(r22, aoi)

r2 <- mosaic(r21, r22, fun = "mean")


# Stack rasters together
rs <- c(r1,r2,r3)
names(rs) <- paste0("r", rep(1:3, each = 4), 1:4)

# create a normalised difference index between all possible combinations
# of two within the raster stack
layers <- names(rs)
for(i in layers){
  print(i)
  layers <- layers[-1]
  for(ii in layers){
    rs[[paste(i, ii, sep = "_")]] <- (rs[[i]]-rs[[ii]])/(rs[[i]]+rs[[ii]])
  }
}


# Transform veg data to spatial points and extract values from raster imagery
pp <- st_as_sf(p2, coords = c("x", "y"), crs = 32634)
d2 <- bind_cols(p2, data.frame(extract(rs, vect(pp))), st_coordinates(pp))
d2 <- d2[,c(7:NCOL(d2))]
d2 <- d2[complete.cases(d2),]

# extract values from raster imagery to the GCG dataset
d3 <- bind_cols(d, data.frame(extract(rs, vect(p))), st_coordinates(p))
d3 <- d3[,c(2,10:NCOL(d3))]
d3 <- d3[complete.cases(d3),]

# extract values from raster imagery to the extra points
d4 <- bind_cols(p3 %>% st_drop_geometry(), data.frame(extract(rs, vect(p3))), st_coordinates(p3))
d4 <- d4[complete.cases(d4),]

# Bind all points together
dall <- bind_rows(d2,d3,d4)
dall <- dall[complete.cases(dall),]


##############################################################
# Fitting the Random Forest classifier

rf <- randomForest(as.factor(cavm) ~ .,
                   data = dall %>% select(-ID, -X, -Y), ntree = 1000)

# OOB estimate of  error rate: 34.5%
# Confusion matrix:
#   b   ds  es   g   w class.error
# b  113   55  60   1   0   0.5065502
# ds  10 2055 464  44   6   0.2031795
# es  20  777 890  51   3   0.4887995
# g    4   97 165 292   8   0.4840989
# w    1   22  16  16 108   0.3374233

# Predict the classes over the whole landscape
rp <- terra::predict(rf, rs, type = "response")

r <- rs[[1]]
r[] <- rp

# Write out the final vegetation classification
writeRaster(r, "output/Mikkuna_vegclass_singlemodel.tif", datatype = "INT1U", overwrite = T)

#Plot the results
plot(r, col = c("burlywood4", "darkolivegreen4", "darkgreen", "chartreuse2", "darkslategray3"))
plot(st_geometry(p), add = T, color = "black", pch = 20)

