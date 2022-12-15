
  # source("~/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/algoritmo_OWM_LiDAR/pintaLAS.R")

  # pkg list
  my.packages <- c("lidR")

  # install and require code
  install.lib <- my.packages[!my.packages %in% installed.packages()]
  for(lib in install.lib) install.packages(lib, dependencies=TRUE)
  sapply(my.packages, require, character=TRUE)
  rm(my.packages, install.lib, lib)

  DirInput <- "/home/felipe/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/datos"

  # DirOutput <- "/home/felipe/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/algoritmo_OWM_LiDAR/resultados"
  DirOutput <- DirInput


# Plot -------------------------------------------------

  # las <- readLAS(file.path(DirInput,"INAER_2011_Alcoy_Core.las"))
  # las <- readLAS(file.path(DirInput,"BABCOCK_2017_Arzua_3B.las"))
  # las <- readLAS(file.path(DirInput,"Guitiriz.las"))
  # las <- readLAS(file.path(DirInput,"Guitiriz_foto.las"))
  # las <- readLAS(file.path(DirInput,"Begonte.las"))
  # las <- readLAS(file.path(DirInput,"V21_group1_densified_point_cloud.las"))
  # las <- readLAS(file.path(DirInput,"VaihingenS09.las"))
  # plot(las)

  # thr <- c(0,2,5,10,15)
  # edg <- c(0, 1.5)
  # chm1 <- grid_canopy(las, 1, pitfree(thr, edg))
  # plot(chm1)


  # semilla <- readLAS(file.path(DirOutput,"INAER_2011_Alcoy_Core_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"BABCOCK_2017_Arzua_3B_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"Guitiriz_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"Guitiriz_foto_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"Begonte_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"V21_group1_densified_point_cloud_salida.las"))
  # semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_min.las"))
  semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_W10_B20_O0.9.las"))
  # semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_W25_B30_O0.9.las"))
  # semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_W28_B30_O0.9.las"))
  # semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_W30_B30_O0.9.las"))
  # semilla <- readLAS(file.path(DirOutput,"VaihingenS09_salida_0.las"))
  plot(semilla)

  # thr <- c(0,2,5,10,15)
  # edg <- c(0, 1.5)
  # chm <- grid_canopy(semilla, 1, pitfree(thr, edg))
  # plot(chm)
