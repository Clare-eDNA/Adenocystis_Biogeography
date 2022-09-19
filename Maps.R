library(maps)
library(mapdata)
library(mapproj)
library(maptools)

maps::map("nzHires", xlim=c(165,179), ylim=c(-50,-35))
nz<-maps::map('nzHires')

samps <- read.csv("2022_Sampling_sites.csv", header=TRUE, stringsAsFactors=T)
samps

palette=(c("#FF0000", "#FF6E00",  "#FFC300", "#FFFF00", "#AAD500", "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082", "#812BA6", "#B857CA", "#D03A87"))
                
par(mar = c(1, 1, 1, 1))
maps::map('nzHires', xlim=c(169,170.5), ylim=c(-46.5,-45))
points(samps$Longitude, samps$Latitude, pch=19, col=palette,cex=1)
pointLabel(samps$Longitude, samps$Latitude, samps$Loctation)

library(ggmap)

map.nz <- map_data(map = "nz")
p1 <- ggplot(map.nz, aes(x = long, y = lat, group=group))
p1 <- p1 + geom_polygon()
p1 <- p1 + labs(title = "New Zealand")+theme_bw()
p1


p2=p1+geom_point(data = samps, aes(x = Longitude, y = Latitude, colour = Location), size = 4,inherit.aes = FALSE, position = "jitter", colour=palette)+theme_bw()+coord_map(xlim=c(169.5, 170.75), ylim=c(-45.5, -46.5))

p2

p2+ geom_text(data = samps, aes(x = Longitude, y = Latitude, label = Location),position = "jitter", hjust = -0.2, colour = "black", size = 3,inherit.aes = FALSE)



#### An interactive map ####
## Shows close-together points a little more easily ##
## Can zoom in 

install.packages(c("leaflet", "sp"))
library(sp)
library(leaflet)
m <- leaflet() %>% setView(lng = 170, lat = -46, zoom = 10)
m %>% addTiles()

samps <- read.csv("2022_Sampling_sites.csv", header=TRUE, stringsAsFactors=T)
samps

# Show first 20 rows from the `quakes` dataset
leaflet(data = samps) %>% addTiles() %>%
  addMarkers(~Longitude, ~Latitude, popup = ~as.character(Loctation), label = ~as.character(Loctation),labelOptions = labelOptions(noHide = T))









## Other useless map code ##




library(ggplot2)
library(ggmap)
library(mapproj)
map <- get_map(location = "New Zealand", zoom = 16)

map.world <- map_data("world")
map.world

par(mar=c(0, 0, 0, 0))
par(mfrow=c(1,1))
world<-map('world')
map(col = 1)



x <- c("ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap")
# install.packages(x) # warning: uncommenting this may take a number of minutes
lapply(x, library, character.only = TRUE) # load the required packages

rm(x)

library(tidyverse)

library(OpenStreetMap)
library(leaflet)


# Part 1 ------------------------------------------------------------------

# Read in the shape file
lnd <- readOGR(dsn = "/Users/clare/Documents/Postdoc/Adenocystis_biogeography/scripts/data/nz-coastlines-and-islands-polygons-topo-150k.shp")

# Look at classes of the data slot
sapply(lnd@data, class)
head(lnd)
tail(lnd)
summary(lnd)
# Change the class of Pop_2001 from dactor to numeric
lnd$Pop_2001 <- as.numeric(lnd$Pop_2001)

plot(lnd)

#####
## OK new package: tmaps ##
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(ggplot2) # tidyverse data visualization package
library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(grid)
library(tmaptools)

# Add fill layer to nz shape
tm_shape(nz) +
  tm_fill() 
## from https://geocompr.robinlovelace.net/adv-map.html#fig:tmshape
map_nz = tm_shape(nz) + tm_polygons()
class(map_nz)
map_nz1 = map_nz +
  tm_shape(nz_elev) + tm_raster(alpha = 0.7)
nz_water = st_union(nz) |> st_buffer(22200) |> 
  st_cast(to = "LINESTRING")
map_nz2 = map_nz1 +
  tm_shape(nz_water) + tm_lines()+tm_compass(type = "8star", position = c("right", "top"))+tm_scale_bar(breaks = c(0, 100, 200), text.size = 1)




nz_region = st_bbox(c(xmin = 1340000, xmax = 1450000,
                      ymin = 5130000, ymax = 5210000),
                    crs = st_crs(nz_height)) |> st_as_sfc()
nz_height_map = tm_shape(nz_elev, bbox = nz_region) +
  tm_raster(style = "cont", palette = "YlGn", legend.show = TRUE) +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 1) +
  tm_scale_bar(position = c("left", "bottom"))
nz_map = tm_shape(nz) + tm_polygons() +
  tm_shape(nz_height) + tm_symbols(shape = 2, col = "red", size = 0.1) + 
  tm_shape(nz_region) + tm_borders(lwd = 3) 
library(grid)
nz_height_map
print(nz_map, vp = viewport(0.8, 0.27, width = 0.5, height = 0.5))



pak_osm <- read_osm(pak, ext=1.1) # reads OSM data based on the bounding box of pak with 10 percent extension

qtm(pak_osm) +
  tm_shape(pak) + 
  tm_fill("Condition", n = 5, palette = "-Blues", colorNA = "grey", 
          textNA = "Missing Values", title = "Access to Water") + 
  tm_borders(alpha = 0.5) + 
  tm_credits("Source: PBS - PSLM 2014 - 15\nDH Corp Ltd.", size = .64, position = c(.62,0.02)) +
  tm_shape(Towers) +
  tm_symbols(color="MNO", shape="MNO", alpha=0.5)



install.packages(c("leaflet", "sp"))
library(sp)
library(leaflet)
m <- leaflet() %>% setView(lng = 170, lat = -46, zoom = 10)
m %>% addTiles()

samps <- read.csv("2022_Sampling_sites.csv", header=TRUE, stringsAsFactors=T)
samps

# Show first 20 rows from the `quakes` dataset
leaflet(data = samps) %>% addTiles() %>%
  addMarkers(~Longitude, ~Latitude, popup = ~as.character(Loctation), label = ~as.character(Loctation),labelOptions = labelOptions(noHide = T))
