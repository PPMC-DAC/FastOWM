# 1.Copyright statement comment -------------------------------------------
# Copyright 2019
# 2.Author comment --------------------------------------------------------
# Sandra Buján
# 3.File description comment ----------------------------------------------
# Function seed points
# 4.source() and library() statements -------------------------------------

  # pkg list
  my.packages <- c("raster","lidR","rgeos","sp","dplyr","ggplot2","tictoc")

  # install and require code
  install.lib <- my.packages[!my.packages %in% installed.packages()]
  for(lib in install.lib) install.packages(lib, dependencies=TRUE)
  sapply(my.packages, require, character=TRUE)
  rm(my.packages, install.lib, lib)

# 5.Parameters --------------------------------------------------
  project.utm <- "+proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs"

  DirInput <- "../datos"
  FileLAS <- "BABCOCK_2017_Arzua_3B.las"
  # FileLAS <- "INSITU_2018_Jenaro.las"
  # FileLAS <- "INAER_2011_Alcoy_Core.las"
  # FileLAS <- "INAER_2011_Alcoy.las"
  # FileLAS <- "sample24.las"

  DirOutput <- "../datos"
  File.Output.LAS <- "BABCOCK_2017_Arzua_3B_seed.las"
  # File.Output.LAS <- "INSITU_2018_Jenaro_seed.las"
  # File.Output.LAS <- "INAER_2011_Alcoy_Core_seed.las"
  # File.Output.LAS <- "INAER_2011_Alcoy_seed.las"
  # File.Output.LAS <- "sample24_seed.las"

  WSize <- 12
  BSize <- 20
  Overlap <- 0.8

# 6.Function definitions --------------------------------------------------
  OWM <- function(PlotPts,WSize,BSize,Overlap){

    Ptos.select <- PlotPts@data[0,]
    Despl <- WSize * (1 - Overlap)
    Despl.Count <- (WSize/Despl) - 1

    Density <- round((PlotPts@header@PHB$`Number of point records`/((PlotPts@bbox[1,2] - PlotPts@bbox[1,1]) * (PlotPts@bbox[2,2] - PlotPts@bbox[2,1]))),2)
    min.Count.Points <- round(0.5 * Density * WSize^2, 0)

    #numero de puntos de la nube
    Npuntos <- PlotPts@header@PHB$`Number of point records`

    ExtentOrg <- extent(PlotPts)
    RasterExtent <- ExtentOrg

    #Selección de mínimos locales con solape entre ventanas
    for (j in c(Despl.Count:0)) {
      for (i in c(Despl.Count:0)) {

        RasterExtent@xmin <- ExtentOrg@xmin - Despl*i
        RasterExtent@xmax <- ExtentOrg@xmax - Despl*i
        RasterExtent@ymin <- ExtentOrg@ymin - Despl*j
        RasterExtent@ymax <- ExtentOrg@ymax - Despl*j

        layout <- raster::raster(RasterExtent, crs = project.utm, resolution = c(WSize,WSize))
        cells <- raster::cellFromXY(layout, PlotPts@data[,c("X","Y")])

        #No seleccionar puntos en aquellas celdas donde no hay puntos suficientes
        cells[cells %in% c(which(sort(table(cells)) < min.Count.Points))] <- NA

        #Seleccionar los puntos con elevación mínima en cada celda
        Ptos.select <- rbind(Ptos.select,PlotPts@data[PlotPts@data[, .I[which.min(Z)], by = cells]$V1,])
        # cat("\nNumero de minimos:", length(PlotPts@data[, .I[which.min(Z)], by = cells]$V1), "\n\n")
        rm(layout,cells)
      }
    }

    # cat("\nLista de mínimos:\n")
    # print(Ptos.select[,1:3])
    # cat("\nNumero total de minimos:", length(Ptos.select$X),'\n\n')

    if (Despl.Count > 0) {
    # Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
    Ptos.select$Score <- 1

    Ptos.select <- aggregate(Ptos.select$Score, as.list(Ptos.select[,]), sum)
    Ptos.select <- Ptos.select[which(Ptos.select$x > 1), 1:ncol(PlotPts@data)]
    # cat("\nCuando hago la eliminacion me quedan",length(Ptos.select$X),"elementos\n\n")
    }

    #Corrección de grandes áreas sin puntos semilla tomando como referencia la dimensión del edificio de mayor tamaño
    if (BSize > 0) {

      layout <- raster::raster(ExtentOrg, crs = project.utm, resolution = c(BSize,BSize))
      cell.out.LLP <- rasterize(Ptos.select[,c("X","Y")], layout,fun="count")
      # cat("Numero de minimos en cada celda:\n",cell.out.LLP[],"\nUn total de",sum(cell.out.LLP[]),"minimos\n\n")
      cells <- raster::cellFromXY(layout, PlotPts@data[,c("X","Y")])
      # cat("¿Qué celda está vacía?\n",c(which(is.na(cell.out.LLP[]))))
      # cat("Numero de celdas vacias:",length(c(which(is.na(cell.out.LLP[])))),'\n')
      cells[cells %in% c(which(!is.na(cell.out.LLP[])))] <- NA
      Ptos.select.add <- PlotPts@data[PlotPts@data[, .I[which.min(Z)], by = cells]$V1,]

      rm(layout,cells)
    }

    # cat("\nAñado",length(Ptos.select.add$X),"minimos\n")
    # print(Ptos.select.add[,1:3])
    # cat('\nEs decir, me he quedado con',length(rbind(Ptos.select, Ptos.select.add)$X))
    # cat(" puntos y he descartado", Npuntos-length(rbind(Ptos.select, Ptos.select.add)$X),"\n\n")

    return(rbind(Ptos.select, Ptos.select.add))
  }

# 7. Run --------------------------------------------------

  LiDAR.Data <- lidR::readLAS(file.path(DirInput,FileLAS))

  # ptm <- proc.time()
  tic()
  Seed.Points <- OWM(LiDAR.Data,  WSize, BSize, Overlap)
  toc()
  # print(proc.time() - ptm)
  Seed.Point.LAS <- LAS(Seed.Points, proj4string = sp::CRS(project.utm), check = TRUE)

  #ggplot(Seed.Points, aes(x=X, y=Y, colour=Classification)) + geom_point() +
  #  scale_colour_continuous(low="forestgreen",high="orange") + coord_fixed(ratio=1) +
  #  ggtitle("LiDAR seed points")

  lidR::writeLAS(Seed.Point.LAS,file.path(DirOutput,File.Output.LAS))


  # las <- readLAS(file.path(DirInput,FileLAS))
  # plot(las)

  # thr <- c(0,2,5,10,15)
  # edg <- c(0, 1.5)
  # chm1 <- grid_canopy(las, 1, pitfree(thr, edg))
  # plot(chm1)

  # semilla <- readLAS(file.path(DirOutput,File.Output.LAS))
  # plot(semilla)
