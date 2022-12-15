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
  # FileLAS <- "INAER_2011_Alcoy.las"
  FileLAS <- "sample24.las"

  DirOutput <- "../datos"
  # File.Output.LAS <- "INAER_2011_Alcoy_seed.las"
  File.Output.LAS <- "sample24_seed.las"

  WSize <- 12
  BSize <- 20
  Overlap <- 0.8

# 6.Function definitions --------------------------------------------------
  OWM <- function(PlotPts,WSize,BSize,Overlap){

    cat('\n')
    print(PlotPts)
    cat("Dimensiones  :",dim(PlotPts@data),'\n')
    cat('\n')
    Ptos.select <- PlotPts@data[0,] #inicializo; vacio
    # print(Ptos.select)
    Despl <- WSize * (1 - Overlap)
    cat("Desplazamiento:   ", Despl, '\n')
    # print(Despl)
    Despl.Count <- (WSize/Despl) - 1
    cat("N desplazamientos:", Despl.Count, '\n')

    Density <- round((PlotPts@header@PHB$`Number of point records`/((PlotPts@bbox[1,2] - PlotPts@bbox[1,1]) * (PlotPts@bbox[2,2] - PlotPts@bbox[2,1]))),2)

    #numero de puntos de la nube
    Npuntos <- PlotPts@header@PHB$`Number of point records`
    cat("\nNumero de puntos: ", Npuntos, "\n\n")
    cat("Print de PlotPts@bbox (bounding-box):\n")
    print(PlotPts@bbox) #tabla 2x2 con los maximos y minimos
    # cat(PlotPts@bbox[1,2],PlotPts@bbox[1,1],PlotPts@bbox[1,2]-PlotPts@bbox[1,1],'\n') #xmax, xmin, xmax-xmin
    # cat(PlotPts@bbox[2,2],PlotPts@bbox[2,1],PlotPts@bbox[2,2]-PlotPts@bbox[2,1],'\n') #ymax, ymin, ymax-ymin
    cat("\nAncho:            ", PlotPts@bbox[1,2]-PlotPts@bbox[1,1],"\n")
    cat("\nAlto:             ", PlotPts@bbox[2,2]-PlotPts@bbox[2,1],"\n")
    cat("\nDensidad:         ", Density,"\n")
    cat("\nEntonces voy a necesitar como mínimo", Npuntos/WSize^2,"celdas\n")

    min.Count.Points <- round(0.5 * Density * WSize^2, 0)
    cat("\nNumero minimo de puntos ", min.Count.Points, "\n\n")

    ExtentOrg <- extent(PlotPts)
    cat("ExtentOrg:\n")
    print(ExtentOrg) #para pintar el extremo minimo y maximo
    cat("\n")
    RasterExtent <- ExtentOrg #esto es simplemente para dotar a RasterExtent de la misma estructura <<extent>>

    #Selección de mínimos locales con solape entre ventanas
    Wcount = 0 #voy a llevar una cuenta de las ventanas
    Ccount = 0 #también la cuenta de todas las celdas
    mcount = 0 #el de minimos
    # plot(ExtentOrg@xmin, ExtentOrg@ymin, col = "white", xlab = "X", ylab = "Y")
    # polygon(x = c(ExtentOrg@xmin,ExtentOrg@xmax,ExtentOrg@xmin,ExtentOrg@xmax),
    #         y = c(ExtentOrg@ymin,ExtentOrg@ymin,ExtentOrg@ymax,ExtentOrg@ymax),
    #         col = "#1b98e0")

    cat("\n///////////////////////////BUCLE///////////////////////////\n")
    for (j in c(Despl.Count:0)) {
      # if(j<0.1) j=0
      for (i in c(Despl.Count:0)) {
        # if(i<0.1) i=0
        # cat(RasterExtent@class) #simplmente para ver que la clase es extent
        RasterExtent@xmin <- ExtentOrg@xmin - Despl*i
        RasterExtent@xmax <- ExtentOrg@xmax - Despl*i
        RasterExtent@ymin <- ExtentOrg@ymin - Despl*j
        RasterExtent@ymax <- ExtentOrg@ymax - Despl*j
        # cat("[",i,",",j,"]\n") #para pintar las coordenadar
        # cat("\n\n\n\n\nVentana" , Wcount <- Wcount+1, "[",i,",",j,"]\n")
        # print(RasterExtent) #para pintar todos los extremos de las diferentes ventanas
        # cat("\n")

        # en RasterExtent tengo los limites min y max de la celda
        layout <- raster::raster(RasterExtent, crs = project.utm, resolution = c(WSize,WSize))
        # raster::raster crea un objeto RasterLayer
        # if(i<1 && j<1)
        # print(layout)
        # tiene esta pinta:
        # class      : RasterLayer
        # dimensions : 6, 10, 60  (nrow, ncol, ncell)
        # resolution : 12, 12  (x, y)
        # extent     : 513748.1, 513868.1, 5403125, 5403197  (xmin, xmax, ymin, ymax)
        # crs        : +proj=utm +zone=32 +ellps=WGS84 +units=m +no_defs


        cells <- raster::cellFromXY(layout, PlotPts@data[,c("X","Y")])
        # raster::cellFromXY es un metodo que obtiene el numero de celdas de un objeto Raster* (en nuestro
        # caso RasterLayer)

        # cat("Posiciones X e Y \n")
        # print(PlotPts@data[,c("X","Y"),])
        # # cat("Dimension: ",dim(c(cells)),'\n')
        # cat('\n')

        # if(i<1 && j<1) {
        #   cat('\n')
        #   print(cells)
        #   cat('\n')
        # }

        #No seleccionar puntos en aquellas celdas donde no hay puntos suficientes;
        #si se cumple la condicion pongo la posicion a 0(NA)

        # print(table(cells))
        # cat('\n')
        # print(sort(table(cells)))
        # cat('\n')
        # cat("Numero de celdas:",length(table(cells)),'\n')
        Ccount <- Ccount + length(table(cells)) #acumulo el numero de celdas

        # Para ilustrar lo que significan los valores de table
        # por ejemplo: compruebo el numero de veces que aparece 59
        # cat(sum(cells==21, na.rm=TRUE))
        # cat("\nOJO:",Npuntos-sum(cells>=1, na.rm=TRUE),"puntos sin celda\n") #para saber el numero de puntos sin celda

        # En cuantas se cumple la condicion?
        # cat("\n¿¿Alguna celda con menos de",min.Count.Points,"puntos??\n")
        # cat(c(which(sort(table(cells)) < min.Count.Points)))
        # cat('\n')
        cells[cells %in% c(which(sort(table(cells)) < min.Count.Points))] <- NA

        # if(i<1 && j<1) {
        #   cat('\n')
        #   print(cells)
        #   cat('\n')
        # }

        #Seleccionar los puntos con elevación mínima en cada celda
        Ptos.select <- rbind(Ptos.select,PlotPts@data[PlotPts@data[, .I[which.min(Z)], by = cells]$V1,])
        #Con .I[] estamos devolviendo la fila en la que se encuentra ese punto;
        #al hacer el PlotPts@data["salida",] me da todas las columnas de cada uno de los minimos
        # cat("\nNumero de fila de los mínimos:\n")
        # print(sort(PlotPts@data[, .I[which.min(Z)], by = cells]$V1)) # row number in PlotPts@data corresponding to each cell
        # Esto me devuelve un data.table con una columna autonombrada V1, por eso es necesario indicarlo!!!
        # cat("\nNumero de minimos:", length(PlotPts@data[, .I[which.min(Z)], by = cells]$V1), "\n\n")
        mcount <- mcount + length(PlotPts@data[, .I[which.min(Z)], by = cells]$V1)
        # Aquí pinto los X Y Z de cada uno de los minimos que he dado; las tres primeras columnas
        # print(PlotPts@data[PlotPts@data[, .I[which.min(Z)], by = cells]$V1,1:3])
        # print(Ptos.select)
        rm(layout,cells)
      }
    }
    cat("\n///////////////////////////FIN BUCLE///////////////////////////\n")

    cat("\nNumero total de celdas :", Ccount,'\n')
    cat("\nNumero total de minimos:", mcount,'\n')
    cat("\nLista de mínimos:\n")
    print(Ptos.select[,1:3])
    cat('\n')


    if (Despl.Count > 0) {
      #Solo lo hago si hay solapamiento; si no lo hay, no puede haber puntos repetidos
      # Únicamente aquellos mínimos locales seleccionados más de una vez se considerarán puntos semilla
      Ptos.select$Score <- 1
      # print(Ptos.select$Score)

      Ptos.select <- aggregate(Ptos.select$Score, as.list(Ptos.select[,]), sum)
      # añade una nueva columna al final del conjunto de datos <<x>>
      # en esta están acumulados los mínimos que son iguales, obteniendo las veces que se obtiene cada uno
      # los puntos iguales se acumulan, ya no tenemos el mismo numero de filas en Ptos.select
      cat("Número de veces que aparece cada mínimo:\n")
      print(Ptos.select$x)
      magrupados <- length(Ptos.select$x)
      cat("\nCuando agrupo los mínimos tengo",magrupados,"elementos\n\n")
      # print(Ptos.select[,ncol(Ptos.select)])
      # print(Ptos.select)
      # print(Ptos.select$X)
      # print(which(Ptos.select$x >= 2))
      Ptos.select <- Ptos.select[c(which(Ptos.select$x >= 2)), 1:ncol(PlotPts@data)]
      # de esta forma me quedo con aquellos que han aparecido mas de una vez
      # y a la vez elimino la última columna
      # print(Ptos.select[,1:3])
      cat("\nCuando hago la eliminacion me quedan",length(Ptos.select$X),"elementos\n\n")
      cat("Es decir, ¡¡he eliminado:",magrupados-length(Ptos.select$X),"puntos porque aparecían solo una vez!!\n\n")
    }

    #Corrección de grandes áreas sin puntos semilla tomando como referencia la dimensión del edificio de mayor tamaño
    if (BSize > 0) {
      # rasterizado de la nube de puntos original
      layout <- raster::raster(ExtentOrg, crs = project.utm, resolution = c(BSize,BSize))
      cell.out.LLP <- rasterize(Ptos.select[,c("X","Y")], layout,fun="count")
      #  rasterize transfer values associated with 'object' type spatial data (points, lines, polygons) to raster cells.
      #  If you want to count the number of points in each grid cell, you can use fun='count'
      # estoy contando cuantos de los minimos seleccionados tengo en las celdas que he obtenido de la nube original
      # print(layout)
      # print(cell.out.LLP)
      cat("Numero de minimos en cada celda:\n",cell.out.LLP[],"\nUn total de",sum(cell.out.LLP[]),"minimos\n\n")
      cells <- raster::cellFromXY(layout, PlotPts@data[,c("X","Y")])
      cat("¿Qué celda está vacía?\n",c(which(is.na(cell.out.LLP[]))))
      if(length(which(is.na(cell.out.LLP[])))==0) cat("Todas las celdas llenas! No tengo que añadir nada\n")
      cat("\n\n")
      cells[cells %in% c(which(is.na(cell.out.LLP[])))] <- NA
      # cells[cells %in% c(which(!is.na(cell.out.LLP[])))] <- NA
      # Es decir, compruebo que celdas están vacías para elegir un mínimo directamente de la nube original PlotPts
      Ptos.select.add <- PlotPts@data[PlotPts@data[, .I[which.min(Z)], by = cells]$V1,]
      # cat("Mínimos adicionales:\n")
      # print(Ptos.select.add[,1:3])
      rm(layout,cells)
    }

    cat("\nPor tanto, inicialmente tengo:",length(Ptos.select$X),"minimos\n")
    cat("\nAñado",length(Ptos.select.add$X),"minimos adicionales\n")
    print(Ptos.select.add[,1:3])
    cat("\nFinalmente tengo",length(rbind(Ptos.select, Ptos.select.add)$X),"minimos\n")
    print(rbind(Ptos.select, Ptos.select.add)[,1:3])
    cat('\n')

    return(rbind(Ptos.select, Ptos.select.add))
  }

# 7. Run --------------------------------------------------

  LiDAR.Data <- lidR::readLAS(file.path(DirInput,FileLAS))

  # elapsed <- system.time(Seed.Points <- OWM(LiDAR.Data,  WSize, BSize, Overlap))
  # print(elapsed)

  tic()
  Seed.Points <- OWM(LiDAR.Data,  WSize, BSize, Overlap)
  toc()

  Seed.Point.LAS <- LAS(Seed.Points, proj4string = sp::CRS(project.utm), check = TRUE)

  # ggplot(Seed.Points, aes(x=X, y=Y, colour=Classification)) + geom_point() +
  #  scale_colour_continuous(low="forestgreen",high="orange") + coord_fixed(ratio=1) +
  #  ggtitle("LiDAR seed points")

  lidR::writeLAS(Seed.Point.LAS,file.path(DirOutput,File.Output.LAS))

# 8. Plot -------------------------------------------------

  # las <- readLAS(file.path(DirInput,FileLAS))
  # plot(las)

  # thr <- c(0,2,5,10,15)
  # edg <- c(0, 1.5)
  # chm1 <- grid_canopy(las, 1, pitfree(thr, edg))
  # plot(chm1)

  # semilla <- readLAS(file.path(DirOutput,File.Output.LAS))
  # plot(semilla)

  # thr <- c(0,2,5,10,15)
  # edg <- c(0, 1.5)
  # chm <- grid_canopy(semilla, 1, pitfree(thr, edg))
  # plot(chm)
