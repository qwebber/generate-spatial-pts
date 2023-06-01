
library(maptools)
library(raster)
library(rgdal)
library(spatstat)
library(rgeos)
library(sp)
library(data.table)
library(ggplot2)

lc<-raster("input/FogoSDSS_RS.tif") # This is a landcover map with different habitat types

## landcover numbers: 
#1 = wetland; 2 = broadleaf; 
#3 = coniferForest; 4 = coniferscrub; 
#5 = mixedwood; 6 = rocky; 
# 7 = water; 8 = lichen; 
# 9 = antrho; 10= N/A

## remove rocky, water, lichen, anthro, and N/A
lc[lc == 10] <- NA
lc[lc == 9] <- NA
lc[lc == 7] <- NA


FogoRoad <- readOGR("input/FogoRoadShapefile/FogoRoads.shp")
FogoRoadShp <- spTransform(FogoRoad, CRS("+proj=utm +zone=21 ellps=GRS80"))
FogoRoadShp <- gBuffer(FogoRoadShp, width = 250)  # create buffer around center

FogoRoadExt <- extract(lc,FogoRoadShp,df = TRUE,cellnumbers = TRUE)

FogoRoadExt <- data.table(FogoRoadExt)[!is.na(FogoSDSS_RS)]
setnames(FogoRoadExt, c("ID","cell","landcover"))

FogoRoadExt[, c("X_coord", "Y_coord") := .(xFromCell(lc, cell),  yFromCell(lc, cell))]
FogoRoadExt[, MasterRowID := .I]

FogoRoadExt
FogoRoadExt[landcover ==  1 | landcover == 2][order(X_coord, Y_coord)]


# Set columns
coord.cols <- c('X_coord', 'Y_coord')
# Create lag and dif column names
lag.cols <- paste('lag', coord.cols, sep = '')
difference.cols <- c('difX', 'difY')
# Use shift  to create lagged cols


landCovFunct <- function(dt, x,y){
  dt2 <- copy(dt)
  dt2 <- dt2[landcover ==  x | landcover == y][order(X_coord, Y_coord)]
  dt2 <- dt2[, shifted := shift(landcover)][landcover != shifted]
  dt2[, rowID := .I][,shiftedID :=shift(rowID)]
  dt2[, (lag.cols) := data.table::shift(.SD, 1, NA, 'lag'),
               .SDcols = coord.cols]
  dt2[, (difference.cols) := .((get(coord.cols[1]) - get(lag.cols[1])) ^2,
                                           (get(coord.cols[2]) - get(lag.cols[2]))^2)]
  dt2[, distLength := sqrt(rowSums(.SD)),
                  .SDcols = difference.cols]
  IDvec <- dt2[distLength < 250 & distLength > 100][sample(1:.N, 6)][,c(rowID,shiftedID)]
  dt2[rowID %in% IDvec, MasterRowID]
}

land_combos <- t(combn(unique(FogoRoadExt$landcover),2))

listOfSites <- lapply(seq_along(land_combos[,1]),FUN = function(i){
        landCovFunct(FogoRoadExt,land_combos[i,1],land_combos[,2])})

aa <- rbindlist(list(listOfSites))

final <- aa[,rbindlist(lapply(.SD, FUN = function(i){
          FogoRoadExt[MasterRowID %in% i]}))]

final[,pairs := rep(1:126,each = 2)]

#holes_SP <- SpatialPointsDataFrame(coords, data = data.frame(holes$holes_ID), proj4string = CRS("+init=epsg:27700"))

coords <- cbind(Easting = as.numeric(as.character(final$X_coord)),
                Northing = as.numeric(as.character(final$Y_coord)))

latLong <- SpatialPoints(coords, proj4string = CRS("+proj=utm +zone=21 ellps=GRS80"))
a1 <- spTransform(latLong, CRS("+proj=longlat +zone=21 ellps=GRS80"))

a1$Long <- coordinates(a1)[, 1]

finalLatLongs <- data.table(coordinates(a1)[, 1], coordinates(a1)[, 2])

dat <- data.table(final$pairs, final$landcoverName)
finalLatLongs <- cbind(finalLatLongs, dat)
colnames(finalLatLongs) <- c("Lat", "Long", "pairs", "landcover")


write.csv(finalLatLongs, "output/finalLatLongs.csv")

Legend<-read.table("input/Legend.csv", header=T, sep=",", quote="",fill=TRUE)
Legend # From this, we can see that different values of our map correspond to different habitats

final[landcover == 1, landcoverName := "Wetland"]
final[landcover == 2, landcoverName := "Broadleaf"]
final[landcover == 3, landcoverName := "ConiferForest"]
final[landcover == 4, landcoverName := "ConiferScrub"]
final[landcover == 5, landcoverName := "MixedWood"]
final[landcover == 6, landcoverName := "Rocky"]
final[landcover == 8, landcoverName := "Lichen"]

shapefile_df <- fortify(FogoRoadShp)

Fogo<-readShapeSpatial("input/FogoMapShapefile/FogoPoly.shp")

legend_title <- "Landcover Types"

# Now the shapefile can be plotted as either a geom_path or a geom_polygon.
# Paths handle clipping better. Polygons can be filled.
# You need the aesthetics long, lat, and group).
png("graphics/Map_snail_sampling.png", width=4000, height=3000, res=500, units="px")
ggplot() +
  scale_fill_manual(legend_title) + 
  geom_polygon(data = Fogo, 
            aes(x = long, y = lat, group = group),fill = 'grey', 
             size = .2) + 
  geom_polygon(data = shapefile_df, 
            aes(x = long, y = lat, group = group),fill = 'white',
            size = .2) +
  geom_path(data = shapefile_df,
            aes(x = long, y = lat, group = group), color = "black",
            size = .2) +
  ylab('Northing') +
  xlab('Easting') +
  geom_point(data = final, 
             aes(x = X_coord, y = Y_coord,color = factor(landcoverName))) + 
  guides(fill = guide_legend(title = NULL)) +
  theme(legend.title = element_blank(),
        legend.key=element_blank(),
        axis.text=element_text(size=7,face = "bold", color = "black"),
        axis.title=element_text(size=7, face = "bold"),
        strip.text = element_text(size=6, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),            
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) #+

dev.off()


