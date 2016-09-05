## Calculos de AIC para las simulaciones

source("~/MEGAsync/bio_comparadaII/bin/R/eee.R")

## Escenario 1
setwd("~/MEGAsync/bio_comparadaII/data/simulations/E1/")
arboles <- c(17, 19, 20, 22, 24, 25, 28, 30, 36, 4, 40, 42, 43, 44, 45,
             49,5,6,9,17,19,20,22,24,25,28,30,36,40,42,43,44,45,49)


for(i in arboles){
  rut_tree <- paste("~/MEGAsync/bio_comparadaII/data/simulations/E1/arbol_", i, sep = "")
  rut_geo <- paste("~/MEGAsync/bio_comparadaII/data/simulations/E1/geo_", i, ".dat", sep = "")
  eeee(rut_tree = rut_tree, rut_geo = rut_geo, max_area = 2)
}
