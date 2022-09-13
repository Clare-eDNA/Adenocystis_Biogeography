library(maps)
library(mapdata)
library(mapproj)
library(maptools)

map('nzHires', xlim=c(165,179), ylim=c(-50,-35))
nz<-map('nzHires')

samps <- read.csv("2022_Sampling_sites.csv", header=TRUE, stringsAsFactors=T)
samps

par(mar = c(1, 1, 1, 1))
map('nzHires', xlim=c(169,170.5), ylim=c(-46.5,-45))
points(samps$Longitude, samps$Latitude, pch=19, col="red",cex=1)
pointLabel(samps$Longitude, samps$Latitude, samps$Loctation)

library(ggmap)

map.nz <- map_data(map = "nz")
p1 <- ggplot(map.nz, aes(x = long, y = lat, group=group))
p1 <- p1 + geom_polygon()
p1 <- p1 + labs(title = "New Zealand")+theme_bw()
p1


p2=p1+geom_point(data = samps, aes(x = Longitude, y = Latitude, colour = Loctation), size = 4,inherit.aes = FALSE, position = "jitter")+theme_bw()+coord_map(xlim=c(169.5, 170.75), ylim=c(-45.5, -46.5))

p2

p2+ geom_text(data = samps, aes(x = Longitude, y = Latitude, label = Loctation),position = "jitter", hjust = -0.2, colour = "black", size = 3,inherit.aes = FALSE)












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
