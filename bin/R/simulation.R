## Simulando arboles bajo diferentes parametros de GEOSSE
library(diversitree)
## Escenarios de simulación
## Escenario I
# especiación extincion igual dispersion igual
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E1")
for(i in 1:50){
  sA <- runif(1, min = 0.1, max = 0.4)
  sB <- sA
  sAB <- 0.5
  xA <- 0.5
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

setwd("~/MEGAsync/bio_comparadaII/data/simulations/E2")
## Escenario II 
for(i in 1:50){
  sA <- runif(1, min = 0.4, max = 0.7)
  sB <- sB
  sAB <- 0.5
  xA <- 0.5
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Recuerde ejecutar sed -i '2d' geo_*
## Para eliminar la segunda fila de los archivos

## Escenario III
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E3")
for(i in 1:50){
  sA <- runif(1, min = 0.8, max = 1)
  sB <- sB
  sAB <- 0.5
  xA <- 0.5
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Recuerde ejecutar sed -i '2d' geo_*
## Para eliminar la segunda fila de los archivos

## Escenario IV
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E4")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- runif(1, min = 0.1, max = 0.4)
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Recuerde ejecutar sed -i '2d' geo_*
## Para eliminar la segunda fila de los archivos

## Escenario V
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E5")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- runif(1, min = 0.4, max = 0.7)
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Escenario VI
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E6")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- runif(1, min = 0.8, max = 1)
  xB <- xA
  dA <- 0.6
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Escenario VII
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E7")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- 0.3
  xB <- 0.3
  dA <- runif(1, min = 0.1, max = 0.4)
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Escenario VIII
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E8")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- 0.3
  xB <- 0.3
  dA <- runif(1, min = 0.4, max = 0.7)
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

## Escenario IX
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E9")
for(i in 1:50){
  sA <- 0.5
  sB <- 0.5
  sAB <- 0.5
  xA <- 0.3
  xB <- 0.3
  dA <- runif(1, min = 0.8, max = 1)
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "01"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "11"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "3 (A B C)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}
