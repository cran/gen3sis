## ----include = FALSE----------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  tidy=T,
  fig.align='center',
  tidy.opts = list(width.cutoff=80),
  results='hold'
)

## ----setup, message=F, eval=F, results='hide'---------------------------------
#  #first we load the required packages
#  library(gen3sis)
#  library(raster)
#  
#  #next we set the working directory to the package path
#  datapath <- system.file("extdata", package="gen3sis")
#  setwd(datapath)

## ----echo=FALSE, message=F, results='hide'------------------------------------
knitr::opts_knit$set(root.dir = '../inst/extdata/')
getwd()
library(gen3sis)
library(raster)

## ----eval=T, message=F--------------------------------------------------------
temperature_brick <- brick('InputRasters/WorldCenter/temp_rasters.grd')
aridity_brick <- brick('InputRasters/WorldCenter/aridity_rasters.grd')
area_brick <- brick('InputRasters/WorldCenter/area_rasters.grd')

## ----eval=TRUE----------------------------------------------------------------
landscapes_list <- list(temp=NULL, arid=NULL, area=NULL)
for(i in 1:nlayers(temperature_brick)){
  landscapes_list$temp <- c(landscapes_list$temp, temperature_brick[[i]])
  landscapes_list$arid <- c(landscapes_list$arid, aridity_brick[[i]])
  landscapes_list$area <- c(landscapes_list$area, area_brick[[i]])
}

## ----eval=TRUE----------------------------------------------------------------
cost_function_null <- function(source, habitable_src, dest, habitable_dest) {
    return(1) #represents cost 4 in a 4degree landscape
}

## ---- echo=FALSE, eval=F------------------------------------------------------
#  landscapes_list_t0 <- list(temp=NULL)
#  landscapes_list_t0$temp <- c(landscapes_list_t0$temp, temperature_brick[[1]])
#  
#  create_input_landscape(landscapes = landscapes_list_t0,
#                                 cost_function = cost_function_null,
#                                 output_directory = file.path(tempdir(), "cost_function_null"),# a directory name to save the files in
#                                 directions = 8, # all surrounding sites from a focus site
#                                 calculate_full_distance_matrices = TRUE  # full distance matrix
#  )
#  

## ---- echo=FALSE, fig.width=7, fig.height=5, fig.cap='This figure shows the connection costs from one site in central Africa to all other sites. To travel to the site in South America that is indicated with the arrow, the travelling cost is 64. The distance matrix was computed using the very simple cost function that has been introduced before and is not adding any penalty.', fig.align='left'----
dist_matrix_null_t0 <- readRDS('CostFunctionExamples/cost_function_null/distances_full/distances_full_0.rds')
landscapes_null <- readRDS('CostFunctionExamples/cost_function_null/landscapes.rds')
landscapes_null_t0 <- na.omit(landscapes_null$temp[, c('x', 'y', '0')])
dist_null_t0_mat <- cbind(landscapes_null_t0[,c('x', 'y')], cost=as.numeric(dist_matrix_null_t0[,140]))
dist_null_t0 <- rasterFromXYZ(dist_null_t0_mat)

maxcost <- 155
mincost <- 0
cost_breaks <- seq(mincost, maxcost, by=1)
cost_colors <- rev(gray(seq(0.03, 0.9, length.out=length(cost_breaks)-1)))

par(mar=c(1,1,1,2))
layout(matrix(c(1,1,1,1,2), ncol=1))

#par(mar=c(1,1,1,2))
image(dist_null_t0, col=cost_colors, breaks=cost_breaks)
#plot(dist_null_t0, col=cost_colors, breaks=cost_breaks, axes=F, box=F, legend.args = list(text = 'connection cost', side = 2, 
#         font = 2, line = 0.5, cex = 0.8))
arrows(10, 8, -54, 8, lwd=3, col='red')
text(-22, 11, labels=paste('connection cost =', round(dist_null_t0_mat$cost[dist_null_t0_mat$x==(-54) & dist_null_t0_mat$y==8], 1)), cex=1.5, font=2)

plot.new()
legend_df <- as.data.frame(cbind(seq(0, length(cost_breaks)-1, length.out=(length(cost_breaks))), rep(0.25, (length(cost_breaks))), cost_breaks))
legend_image <- rasterFromXYZ(legend_df, res=0.01)
plot(legend_image, legend.only=T, col=cost_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.45, 0.6), 
   axis.args=list(at=seq(mincost, maxcost, 25),labels=seq(mincost, maxcost, 25)), legend.args=list(text='connection cost'))

## ----eval=TRUE----------------------------------------------------------------
cost_function_water <- function(source, habitable_src, dest, habitable_dest) {
  if(!all(habitable_src, habitable_dest)) {
    return(2)
  } else {
    return(1)
  }
}

## ---- echo=FALSE, eval=F------------------------------------------------------
#  create_input_landscape(landscapes = landscapes_list_t0,
#                                 cost_function = cost_function_water,
#                                 output_directory = file.path(tempdir(), "cost_function_water"),# a directory name to save the files in
#                                 directions = 8, # all surrounding sites from a focus site
#                                 calculate_full_distance_matrices = TRUE  # full distance matrix
#  )
#  

## ---- echo=FALSE, fig.width=7, fig.height=5, fig.align='left'-----------------
dist_matrix_water_t0 <- readRDS('CostFunctionExamples/cost_function_water/distances_full/distances_full_0.rds')
landscapes_water <- readRDS('CostFunctionExamples/cost_function_water/landscapes.rds')
landscapes_water_t0 <- na.omit(landscapes_water$temp[, c('x', 'y', '0')])
dist_water_t0_mat <- cbind(landscapes_water_t0[,c('x', 'y')], cost=as.numeric(dist_matrix_water_t0[,140]))
dist_water_t0 <- rasterFromXYZ(dist_water_t0_mat)

cmaxcost <- 155
mincost <- 0
cost_breaks <- seq(mincost, maxcost, by=1)
cost_colors <- rev(gray(seq(0.03, 0.9, length.out=length(cost_breaks)-1)))

par(mar=c(1,1,1,2))
layout(matrix(c(1,1,1,1,2), ncol=1))

#par(mar=c(1,1,1,2))
image(dist_water_t0, col=cost_colors, breaks=cost_breaks)
#plot(dist_water_t0, col=cost_colors, breaks=cost_breaks, axes=F, box=F, legend.args = list(text = 'connection cost', side = 2, 
#         font = 2, line = 0.5, cex = 0.8))
arrows(10, 8, -54, 8, lwd=3, col='red')
text(-22, 11, labels=paste('connection cost =', round(dist_water_t0_mat$cost[dist_water_t0_mat$x==(-54) & dist_water_t0_mat$y==8], 1)), cex=1.5, font=2)

plot.new()
legend_df <- as.data.frame(cbind(seq(0, length(cost_breaks)-1, length.out=(length(cost_breaks))), rep(0.25, (length(cost_breaks))), cost_breaks))
legend_image <- rasterFromXYZ(legend_df, res=0.01)
plot(legend_image, legend.only=T, col=cost_colors, horizontal=T, smallplot=c(0.2, 0.8, 0.45, 0.6), 
   axis.args=list(at=seq(mincost, maxcost, 25),labels=seq(mincost, maxcost, 25)), legend.args=list(text='connection cost'))

## ----eval=F-------------------------------------------------------------------
#  create_input_landscape(landscapes = landscapes_list,
#               cost_function = cost_function_water,
#               output_directory = file.path(tempdir(), "WorldCenter"),
#               timesteps = paste0(round(seq(150, 100, length.out = 301),2), "Ma"),
#               calculate_full_distance_matrices = F)

