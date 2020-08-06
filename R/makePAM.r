#' Presence-absence matrices
#'
#' Creates a presence-absence matrix (PAM) from an sf object
#' @param sf_object an object of class "sf" containing all of the polygon data to be converted to a PAM
#' @param taxon_names column name in the sf object specifying taxon names to be used to create the PAM
#' @param resolution required resolution for the PAM. The appropriate size depends on your data and question. Default is 100 (km) for a Berhman equal area grid (see `crs` below).
#' @param sampling one of "conservative" (the default but not necessarily the best), "liberal", or a numeric value from 1:100
#' @param file_name path to output file to write the PAM to file as R data (.rds). If no path is specified, the file will not be written.
#' @param crs coordinate reference system for the PAM. Default is "+proj=cea +lat_ts=30 +units=km" (a Berhman equal area projection)
#' @param clip_to_world One of "none", "sampling", or "liberal". None does no clipping and returns the full PAM. "sampling" and "liberal" sets grid cells with no land to NA.
#' @param returnPAM logical. If `TRUE` the PAM is returned and can be used direction.
#' @details TBC
#' @return A list of two or three elements (if `clip_to_world=TRUE`). The first element is the clipped PAM (if `clip_to_world=TRUE`) or the unclipped PAM (if `clip_to_world=FALSE`) and the last element is an empty raster object with the same crs and extent as the PAM. This is needed for plotting the PAM as e.g. a species richness map.
#' @examples  
#'
#' @export
sf_to_pam <- function (sf_object, taxon_names=NULL, resolution=100, sampling="conservative", file_name=NULL, crs="+proj=cea +lat_ts=30 +units=km", clip_to_world="none", returnPAM=TRUE) {

require(sf)
require(dplyr)
require(fasterize)
require(raster)
require(maptools)
require(velox)
require(stars)
require(maptools)

    if (is.null(file) & returnPAM==FALSE) { stop("You have set file=NULL and returnPAM=FALSE you daft apeth. You must define an output file, set returnPAM=TRUE, or both if you want to output a PAM.")}

    if (is.null(taxon_names)) {
        stop("Column of taxon_names must be defined as a character in the function call.")
    }

   print("Projecting shape data...")
    maps_sf <- st_transform(sf_object, crs=crs)
    empty_behrman_raster <- raster(maps_sf, res=resolution)
    crs(empty_behrman_raster) <- crs

    
        maps_sf$count <- 1
        grid <- st_as_stars(st_bbox(empty_behrman_raster), values=0, nx=empty_behrman_raster@ncols, ny=empty_behrman_raster@nrows,
            xlim=c(xmin(empty_behrman_raster), xmax(empty_behrman_raster)), ylim=c(ymin(empty_behrman_raster), ymax(empty_behrman_raster)))
        

    if (is.numeric(sampling)) {
       null_rast <- disaggregate(empty_behrman_raster, fact=c(10,10), method="")
       crs(null_rast) <- crs
        }

    map_names <- unique(as.data.frame(maps_sf)[,taxon_names])
    n_species <- length(map_names)

    pres_ab <- matrix(NA, ncol=n_species, nrow=ncell(empty_behrman_raster))
    colnames(pres_ab) <- map_names

    print("Projection complete. Creating map for each taxon.")
    for (i in 1:n_species) {

        current_map <- maps_sf %>% filter(.data[[taxon_names]]==map_names[i])   #####

        if (dim(current_map)[1]>1) {
            current_map <- st_cast(current_map, "MULTIPOLYGON")
            }

        if (sampling=="conservative") {
            sp_rast <- fasterize(current_map, empty_behrman_raster, field=NULL)
            idx <- !is.na(sp_rast@data@values)
            pres_ab[idx, i] <- rep(1, sum(idx))
            }

        if (sampling=="liberal") {
            current_map <- st_rasterize(current_map["count"], grid, options = c("ALL_TOUCHED=TRUE", "MERGE_ALG=REPLACE"))
            sp_rast <- fasterize(st_as_sf(current_map), raster=empty_behrman_raster, field="count")
            idx <- !is.na(sp_rast@data@values) & sp_rast@data@values > 0
            pres_ab[idx, i] <- rep(1, sum(idx))
            }

        if (is.numeric(sampling)) {
            sp_rast <- fasterize(current_map, null_rast, field=NULL, background=0)
            vx <- velox(sp_rast)
            vx$aggregate(factor=c(10,10), aggtype="sum")
            sp_rast <- vx$as.RasterLayer(band=1)
            
            idx <- which(sp_rast@data@values >= sampling)
            pres_ab[idx, i] <- rep(1, length(idx))
        }

        cat(i, "of ", n_species, "\r") 
        flush.console()
  }

    if (sampling=="conservative" & clip_to_world=="sampling") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        worldrast<-fasterize(world, empty_behrman_raster, field=NULL)

        pres_ab_clip <- pres_ab
        pres_ab_clip[is.na(worldrast@data@values)] <- NA 
        }

    if (sampling=="liberal" & clip_to_world=="sampling") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        world$count <- 1
        #grid <- st_as_stars(st_bbox(world))

        worldrast <- st_rasterize(world["count"], grid, options = c("ALL_TOUCHED=TRUE", "MERGE_ALG=REPLACE"))
        worldrast<-fasterize(st_as_sf(worldrast), empty_behrman_raster, field="count")

        pres_ab_clip <- pres_ab
        pres_ab_clip[worldrast@data@values==0] <- NA   
        }
    
    if (is.numeric(sampling) & clip_to_world=="sampling") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        worldrast<-fasterize(world, null_rast, field=NULL)

        worldrast <- aggregate(worldrast, fact=10, sum)
        worldrast@data@values[worldrast@data@values<sampling] <- NA

        pres_ab_clip <- pres_ab
        pres_ab_clip[is.na(worldrast@data@values)] <- NA 
        }

    if (sampling=="liberal" & clip_to_world=="liberal") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        world$count <- 1

        worldrast <- st_rasterize(world["count"], grid, options = c("ALL_TOUCHED=TRUE", "MERGE_ALG=REPLACE"))
        worldrast<-fasterize(st_as_sf(worldrast), empty_behrman_raster, field="count")

        pres_ab_clip <- pres_ab
        pres_ab_clip[worldrast@data@values==0] <- NA   
        }

    if (sampling=="conservative" & clip_to_world=="liberal") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        world$count <- 1

        worldrast <- st_rasterize(world["count"], grid, options = c("ALL_TOUCHED=TRUE", "MERGE_ALG=REPLACE"))
        worldrast<-fasterize(st_as_sf(worldrast), empty_behrman_raster, field="count")

        pres_ab_clip <- pres_ab
        pres_ab_clip[worldrast@data@values==0] <- NA   
        }

    if (is.numeric(sampling) & clip_to_world=="liberal") {
        data("wrld_simpl", package = 'maptools')
        world <- st_as_sf(wrld_simpl)
        world <- st_transform(world, crs=crs)
        world$count <- 1

        worldrast <- st_rasterize(world["count"], grid, options = c("ALL_TOUCHED=TRUE", "MERGE_ALG=REPLACE"))
        worldrast<-fasterize(st_as_sf(worldrast), empty_behrman_raster, field="count")

        pres_ab_clip <- pres_ab
        pres_ab_clip[worldrast@data@values==0] <- NA   
        }

    if (!is.null(file_name)) {
        print(paste("PAM complete, writing output to: ", file_name, ".rds", sep=""))
        saveRDS(list(pres_ab, empty_behrman_raster) , file=paste(file_name,  ".rds", sep=""))

        if (clip_to_world=="sampling" || clip_to_world=="liberal"){
            print(paste("PAM complete, writing output to: ", file_name, "_clipped.rds", sep="")  )
            saveRDS(list(pres_ab, pres_ab_clip, empty_behrman_raster) , file=paste(file_name,  "_clipped.rds", sep=""))
            }
        }

  if (returnPAM) {
        if(clip_to_world=="none") {return(list(pres_ab, empty_behrman_raster))}
        if(clip_to_world=="sampling" || clip_to_world=="liberal") {return(list(pres_ab, pres_ab_clip, empty_behrman_raster))}
    }

}





#' Make raster object
#'
#' Creates a raster from a PAM
#' @param pam a pam object produced by sf_to_pam
#' @param plot_map logical. If `TRUE` plot the raster as a map
#' @param clipped logical. If `TRUE` will rasterise/plot a clipped PAM if a clipped PAM is present.
#' @details TBC
#' @return A a raster object.
#' @examples 
#'
#' @export
pam_to_raster <- function(pam, plot_map=TRUE, clipped=TRUE) {
  
require(raster)
require(rasterVis)
require(ggplot2)
require(ggspatial)
  
  if (length(pam)==3 & clipped==TRUE) { pam_rast <- pam[[2]]  }
  if (length(pam)==3 & clipped==FALSE) { pam_rast <- pam[[1]]  }
  if (length(pam)==2) { pam_rast <- pam[[1]]  }
  
PAM_raster <- pam[[length(pam)]]
values(PAM_raster) <- rowSums(pam_rast, na.rm=TRUE)
PAM_raster@data@values[PAM_raster@data@values==0] <- NA 

if (plot_map==TRUE) {
  p <- ggplot() +  
  layer_spatial(PAM_raster) +
  scale_fill_continuous(name="Value", type = "viridis", na.value = NA) +
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 20))
  plot(p)
}

return(PAM_raster)
}
