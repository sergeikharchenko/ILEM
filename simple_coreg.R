setwd("C:/ad") # to edit
list.files()

necessary_packages <- c("future", "sp", "terra", "raster", "rgeos", "furrr", "Rsagacmd")
install.packages(necessary_packages[!(necessary_packages %in% installed.packages()[,1])])


library(future)
library(sp)
library(terra)
library(raster)
library(rgeos)
library(furrr)
library(Rsagacmd)

# PARALLEL PROCESSING SETTINGS
options(future.globals.maxSize= 3048^3)
cl <- makeClusterPSOCK(availableCores() - 6) # you should consider free RAM for each CPU process
plan(cluster, workers = cl)

# DEM FILES NAMES
ref_name <- "elev20120929.tif" # to edit
sou_name <- "elev20220714.tif" # to edit

# READ THE DATA
(r_ref <- rast(ref_name))
(r_sou <- rast(sou_name))
#r_ref <- resample(r_ref, rast(ext = ext(r_ref), res = 0.3, crs = crs(r_ref)))
#r_sou <- resample(r_sou, rast(ext = ext(r_sou), res = 0.3, crs = crs(r_sou)))
#r_sou <- crop(r_sou, r_ref)

# MAKE THE NETWORK
size <- 300 # to edit

br1 <- boundaries(resample(r_ref, rast(ext = ext(r_ref), res = size/10)), inner=T)
br1[br1 == 0] <- 1
br1 <- as.polygons(br1)

tryCatch({
  ss1 <- spsample(x = as(br1, "Spatial"), cellsize = size, type = "hexagonal")
}, error = function(e) {assign(ss1, spsample(x = as(as.polygons(ext(br1)), "Spatial"), cellsize = size, type = "hexagonal"))})

writeVector(vect(ss1), "ss.shp", overwrite = TRUE)

br2 <- boundaries(resample(r_sou, rast(ext = ext(r_sou), res = size/10)), inner=T)
br2[br2 == 0] <- 1
br2 <- as.polygons(br2)

tryCatch({
  ss1 <- ss1[relate(vect(ss1), br2, "intersects")]
}, error = function(e) {assign(ss1, ss1[relate(vect(ss1), as.polygons(ext(br2)), "intersects")] )})

hex_grid <- HexPoints2SpatialPolygons(ss1, dx = 2*size) # 2 - to edit
hex_grid$dx <- rep(0,length(hex_grid))
hex_grid$dy <- rep(0,length(hex_grid))
hex_grid$err <- rep(0,length(hex_grid))
hex_grid_v <- vect(hex_grid)
writeVector(hex_grid_v, "hex_grid_v.shp", overwrite = T)

# PLOT THE NETWORK
par(mfrow = c(1,1))
plot(r_ref)
crs(ss1) <- crs(r_ref)
points(ss1)
plot(hex_grid, add = T)

# SHIFT STEPS FOR EACH ITERATIVE
seq_list <- list(seq(-15, 15, 5), seq(-8, 8, 2), seq(-3.2, 3.2, 0.4)) # to edit

# COMPUTATION FOR COARSE SHIFT PARAMETERS
time1 <- Sys.time()
.r_ref <- wrap(r_ref)
.r_sou <- wrap(r_sou)
.hex_grid_v <- wrap(hex_grid_v)

command_set <- vector(mode = "list", length = length(hex_grid))
for (i in 1:length(hex_grid)) {
  if (i %% length(cl) == 0) print(round(100 * i / length(hex_grid), 1))
  command_set[[i]] <- future({
    hgv <- vect(.hex_grid_v)[i,]
    r_ref_clip <- mask(crop(rast(.r_ref), hgv), hgv)
    r_sou_clip <- mask(crop(rast(.r_sou), hgv), hgv)
    
    expg1 <- expand.grid(hex_grid@data[i,1] + seq_list[[1]], hex_grid@data[i,2] + seq_list[[1]])
    plan_err <- rep(NA, nrow(expg1))
    for (j in 1:nrow(expg1)) {
      plan_err[j] <- global(resample(shift(r_sou_clip, expg1[j,1], expg1[j,2]), r_ref_clip) - r_ref_clip, function(x) sd(x, na.rm = T))[[1]]
    }
    c(unlist(expg1[which.min(plan_err),]), min(plan_err, na.rm = T))
  })
}
output <- value(command_set)
output <- as.data.frame(do.call(rbind, output))
hex_grid@data <- output
colnames(hex_grid@data) <- c("dx","dy", "err")

# PLOT THE FIRST ITERATIVE RESULTS
par(mfrow = c(2,2))
plot(density(hex_grid$err, na.rm = T), main = "sd_density")
hist(hex_grid$err, main = "sd_hist")
plot(density((hex_grid$dx^2 + hex_grid$dy^2)^0.5), main = "shift_density")
hist((hex_grid$dx^2 + hex_grid$dy^2)^0.5, main = "shift_hist")
par(mfrow = c(1,1))

# FILTERING NETWORK
ind1 <- (hex_grid$dx^2 + hex_grid$dy^2)^0.5 < 6 & hex_grid$err < 1.5 # to edit
ind1 <- ifelse(is.na(ind1), F, ind1 )
plot(r_ref)
points(ss1, col = ind1)

ss1 <- ss1[ind1,]
hex_grid_v <- hex_grid_v[ind1,]
.hex_grid_v <- wrap(hex_grid_v)
hex_grid <- hex_grid[ind1,]

# THE SECOND AND FURTHER 

for (sh in 2:length(seq_list)) {
  #.r_ref <- wrap(r_ref)
  #.r_sou <- wrap(r_sou)
  
  command_set <- vector(mode = "list", length = length(hex_grid))
  for (i in 1:length(hex_grid)) {
    if (i %% length(cl) == 0) print(round(100 * i / length(hex_grid), 1))
    command_set[[i]] <- future({
      hgv <- vect(.hex_grid_v)[i,]
      r_ref_clip <- mask(crop(rast(.r_ref), hgv), hgv)
      r_sou_clip <- mask(crop(rast(.r_sou), hgv), hgv)
      
      expg1 <- expand.grid(hex_grid@data[i,1] + seq_list[[sh]], hex_grid@data[i,2] + seq_list[[sh]])
      plan_err <- rep(NA, nrow(expg1))
      for (j in 1:nrow(expg1)) {
        plan_err[j] <- global(resample(shift(r_sou_clip, expg1[j,1], expg1[j,2]), r_ref_clip) - r_ref_clip, function(x) sd(x, na.rm = T))[[1]]
      }
      c(unlist(expg1[which.min(plan_err),]), min(plan_err, na.rm = T))
    })
  }
  output <- value(command_set)
  hex_grid@data <- as.data.frame(do.call(rbind, output))
  newcoords <- coordinates(hex_grid) + hex_grid@data[,1:2]
}
time2 <- Sys.time()
time2 - time1

crs(ss1) <- crs(hex_grid) <- crs(hex_grid_v) <- crs(r_ref)

writeVector(vect(ss1), paste0("ss",size,"orig.shp"), overwrite=TRUE)
writeVector(vect(SpatialPoints(coordinates(hex_grid) - hex_grid@data[,1:2])), paste0("ss",size,"sh.shp"), overwrite=TRUE)
writeVector(hex_grid_v, paste0("hex_grid_",size,"orig.shp"), overwrite=TRUE)

orig_points <- vect(paste0("ss",size,"orig.shp"))
crs(orig_points) <- crs(r_ref)
sh_points <- vect(paste0("ss",size,"sh.shp"))
crs(sh_points) <- crs(r_ref)
hg <- vect(paste0("hex_grid_",size,"orig.shp"))
crs(hg) <- crs(r_ref)

plot(orig_points)
points(sh_points, col = "red", pch = 19)
plot(hg, add = T)

writeVector(orig_points, "ss_all_orig.shp", overwrite = T)
writeVector(sh_points, "ss_all_sh.shp", overwrite = T)

parallel::stopCluster(cl)

sagaversion <- "9.3.1" # to edit, saga version

if (!paste0("saga-",sagaversion,"_x64") %in% list.files()) {
  download.file(paste0("https://sourceforge.net/projects/saga-gis/files/SAGA%20-%20",substr(sagaversion, 1, 1),"/SAGA%20-%20",sagaversion,"/saga-",sagaversion,"_x64.zip"), destfile = "saga_gis.zip", method = "libcurl")
  unzip("saga_gis.zip")
  unlink("saga_gis.zip")
}

saga <- saga_gis(saga_bin = paste0("./saga-",sagaversion,"_x64/saga_cmd.exe"),
                 raster_format = "GeoTIFF", 
                 vector_format = "ESRI Shapefile")
system(paste0('./saga-',sagaversion,'_x64/saga_cmd.exe pj_georeference 1 -REF_SOURCE="ss_all_sh.shp" -REF_TARGET="ss_all_orig.shp" -GRID="',sou_name,'" -TARGET_OUT_GRID="',paste0(strsplit(sou_name, "\\.")[[1]][1],"_ref.tif"),'"'))
#saga$io_gdal$export_geotiff(grids = paste0(strsplit(sou_name, "\\.")[[1]][1],"_ref.sdat"), 
#                            file = paste0(strsplit(sou_name, "\\.")[[1]][1],"_ref.tif"))


## DEM COREG
(h_err <- unlist(sapply(1:length(hg), function(x) global(crop(r_ref, hg[x,]) - resample(crop(r_sou, hg[x,]), crop(r_ref, hg[x,])), function(a) median(a, na.rm = T)))))
sd(h_err, na.rm = T)
writeRaster(resample(r_sou, r_ref) - r_ref, "rdiff_initial.tif", overwrite=TRUE)


r_sou <- rast(paste0(strsplit(sou_name, "\\.")[[1]][1],"_ref.tif"))
(h_err <- unlist(sapply(1:length(hg), function(x) global(crop(r_ref, hg[x,]) - resample(crop(r_sou, hg[x,]), crop(r_ref, hg[x,])), function(a) median(a, na.rm = T)))))
sd(h_err, na.rm = T)
writeRaster(resample(r_sou, r_ref) - r_ref, "rdiff_planform.tif", overwrite=TRUE)


#
r_error <- rast(kriging::kriging(x = crds(orig_points)[,1], y = crds(orig_points)[,2], response = h_err, lags = 2)$map, type="xyz")

#dferrors <- data.frame(x = crds(orig_points)[,1], y = crds(orig_points)[,2], response = h_err)
#lmerrors <- lm(response ~ ., dferrors)
#rerrors <- rast(ext = ext(r_ref), res = 3*res(r_ref), crs = crs(r_ref))
#rerrors[] <- 1
#newdata <- as.data.frame(crds(as.points(rerrors)))
#r_error <- rast(cbind(newdata, predict(lmerrors, newdata)), type="xyz")

plot(r_error)
writeRaster(r_error, "r_error.tif", overwrite=TRUE)


crs(r_error) <- crs(r_sou)
r_sou <- r_sou + resample(r_error, r_sou)
writeRaster(r_sou, paste0(r_sou@ptr@.xData$names,"_zcor.tif"), overwrite=TRUE)

(rand_err <- sd(unlist(sapply(1:length(hg), function(x) global(crop(r_ref, hg[x,]) - resample(crop(r_sou, hg[x,]), crop(r_ref, hg[x,])), function(a) median(a, na.rm = T))))))

r_diff <- resample(r_sou, r_ref) - r_ref
r_diff <- trim(r_diff)
plot(r_diff)
writeRaster(r_diff, "r_diff.tif", overwrite=TRUE)

sum(r_diff[r_diff > 0], na.rm = T) * prod(res(r_diff))
sum(r_diff[r_diff < 0], na.rm = T) * prod(res(r_diff))

sum(r_diff[r_diff > 3*abs(rand_err)], na.rm = T) * prod(res(r_diff))
sum(r_diff[r_diff < 3*-abs(rand_err)], na.rm = T) * prod(res(r_diff))