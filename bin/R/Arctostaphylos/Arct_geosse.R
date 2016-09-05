## Calculo AIC datos empiricos

library(diversitree)
library(BioGeoBEARS)
arbol <- read.tree("~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/output/Arctostaphylos.new") 
geogfn = "~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/Arctostaphylos_geo.dat"

datos <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
plot(arbol)
arbol$tip.label
datos <- as.data.frame(datos@df)
datos1 <- vector()
for( i in 1:nrow(datos)){
  if(datos$C[i] ==1 & datos$F[i] ==0){
    datos1[i] <-1
  }
  if(datos$C[i] ==1 & datos$F[i] ==1 ){
    datos1[i] <-0
  }
  if(datos$F[i] ==1 & datos$C[i] == 0){
    datos1[i] <-2
  }
}
names(datos1) <- rownames(datos)
length(datos1)
lik1 <- make.geosse(arbol, states =  datos1)
pars <- starting.point.geosse(arbol)
likelihod <- find.mle(lik1, x.init = pars )

likelihod$lnLik

## Calculo valor AIC
aic_muse_arct <- -2 * (likelihod$lnLik) + 2 * 6
