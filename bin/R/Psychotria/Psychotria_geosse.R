## Calculo AIC datos empiricos

library(diversitree)
library(BioGeoBEARS)
arbol <- read.tree("~/MEGAsync/bio_comparadaII/result/Psychotria/final/output/Psychotria.new") 
geogfn = "~/MEGAsync/bio_comparadaII/result/Psychotria/final/Psychotria_geo.dat"

datos <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
plot(arbol)
arbol$tip.label
datos <- as.data.frame(datos@df)
datos1 <- vector()
for( i in 1:nrow(datos)){
  if(datos$K[i] ==1){
    datos1[i] <-1
  }
  if(datos$O[i] ==1){
    datos1[i] <-2
  }
  if(datos$M[i] ==1){
    datos1[i] <-3
  }
  if(datos$H[i] ==1){
    datos1[i] <-4
  }
}
names(datos1) <- rownames(datos)
length(datos1)
lik1 <- make.musse(arbol, states =  datos1, k = 4)
pars <- starting.point.musse(arbol, k = 4)

likelihod <- find.mle(lik1, x.init = pars )
str(likelihod)
likelihod$lnLik

## Calculo valor AIC
aic_muse_psyco <- -2 * (likelihod$lnLik) + 2 * 20
