library(maps)
library(mapdata)
library(mapproj)
library(maptools)

map('nzHires', xlim=c(165,179), ylim=c(-50,-35))
nz<-map('nzHires')

samps <- read.csv("2022_Sampling_sites.csv", header=TRUE, stringsAsFactors=T)
samps
