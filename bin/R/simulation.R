## Simulando arboles bajo diferentes parametros de GEOSSE
library(diversitree)
## Escenarios de simulación
## Escenario I
# Parametros Simetricos
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E1")
for(i in 1:50){
  sA <- 2.0
  sB <- sA
  sAB <- 0
  xA <- 1.9
  xB <- xA
  dA <- 2.0
  dB <- dA
  pars <- c(sA, sB, sAB, xA, xB, dA, dB)
  arbol <- tree.geosse(pars=pars, max.taxa = 40)
  if(is.null(arbol)==FALSE){
    print(i)
    write.tree(arbol, file = paste("arbol_", i, sep = ""))
    p <- starting.point.geosse(arbol)
    lik1 <- make.geosse(tree =arbol, states = arbol$tip.state)
    ml1 <- find.mle(lik1, x.init = p)
    salida <- data.frame()
    salida[1,"LnL"] <- ml1$lnLik
    salida[1, "AIC"]<- -2 * (ml1$lnLik) + 2 * 5
    write.csv(salida, file = paste("geosse_", i, sep = ""))
    datos <- data.frame()
    for(a in 1:length(arbol$tip.label)){
      datos[a,1] <- arbol$tip.label[a]
      if(arbol$tip.state[a]==0){
        datos[a,2] <- "11"
      }
      if(arbol$tip.state[a]==1){
        datos[a,2] <- "10"
      }
      if(arbol$tip.state[a]==2){
        datos[a,2] <- "01"
      }
    }
    archivo <- file(paste("geo_", i, ".dat", sep = ""))
    writeLines(paste(length(arbol$tip.label), "2 (A B)", sep = "\t"), archivo)
    close(archivo)
    write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
  }
}

# 
# setwd("~/MEGAsync/bio_comparadaII/data/simulations/E2")
# ## Escenario II 
# ## Especiacion asimetrica
# for(i in 1:50){
#   sA <- 0.5
#   sB <- 2.0
#   sAB <- 0
#   xA <- 1.9
#   xB <- 1.9
#   dA <- 2.0
#   dB <- 2.0
#   pars <- c(sA, sB, sAB, xA, xB, dA, dB)
#   arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
#   if(is.null(arbol)==FALSE){
#     print(i)
#     write.tree(arbol, file = paste("arbol_", i, sep = ""))
#     datos <- data.frame()
#     for(a in 1:length(arbol$tip.label)){
#       datos[a,1] <- arbol$tip.label[a]
#       if(arbol$tip.state[a]==0){
#         datos[a,2] <- "10"
#       }
#       if(arbol$tip.state[a]==1){
#         datos[a,2] <- "01"
#       }
#       if(arbol$tip.state[a]==2){
#         datos[a,2] <- "11"
#       }
#     }
#     archivo <- file(paste("geo_", i, ".dat", sep = ""))
#     writeLines(paste(length(arbol$tip.label), "2 (A B)", sep = "\t"), archivo)
#     close(archivo)
#     write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
#   }
# }
# 
# ## Recuerde ejecutar sed -i '2d' geo_*
# ## Para eliminar la segunda fila de los archivos
# 
# ## Escenario III
# ## Extincion asimetrica
# setwd("~/MEGAsync/bio_comparadaII/data/simulations/E3")
# for(i in 1:50){
#   sA <- 2.0
#   sB <- 2.0
#   sAB <- 0
#   xA <- 1.9
#   xB <- 0.4
#   dA <- 2.0
#   dB <- 2.0
#   pars <- c(sA, sB, sAB, xA, xB, dA, dB)
#   arbol <- tree.geosse(pars=pars, max.taxa = as.integer(runif(1, min = 30, max = 70)))
#   if(is.null(arbol)==FALSE){
#     print(i)
#     write.tree(arbol, file = paste("arbol_", i, sep = ""))
#     datos <- data.frame()
#     for(a in 1:length(arbol$tip.label)){
#       datos[a,1] <- arbol$tip.label[a]
#       if(arbol$tip.state[a]==0){
#         datos[a,2] <- "10"
#       }
#       if(arbol$tip.state[a]==1){
#         datos[a,2] <- "01"
#       }
#       if(arbol$tip.state[a]==2){
#         datos[a,2] <- "11"
#       }
#     }
#     archivo <- file(paste("geo_", i, ".dat", sep = ""))
#     writeLines(paste(length(arbol$tip.label), "2 (A B)", sep = "\t"), archivo)
#     close(archivo)
#     write.table(datos, file = paste("geo_", i, ".dat", sep = ""), row.names = F, quote = FALSE, sep = "\t", append = TRUE)
#   }
# }
# 
# ## Recuerde ejecutar sed -i '2d' geo_*
# ## Para eliminar la segunda fila de los archivos
