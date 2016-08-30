## Covirtiendo arboles a formato newick

library(ape)
arbol <- read.nexus("~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/output/Arctostaphylos.tree")
write.tree(arbol, "~/MEGAsync/bio_comparadaII/result/Arctostaphylos/final/output/Arctostaphylos.new")

arbol <- read.nexus("~/MEGAsync/bio_comparadaII/result/Ceanothus/final/output/Ceanothus.tree")
write.tree(arbol, "~/MEGAsync/bio_comparadaII/result/Ceanothus/final/output/Ceanothus.new")

arbol <- read.nexus("~/MEGAsync/bio_comparadaII/result/Psychotria/final/output/Psychotria.tree")
write.tree(arbol, "~/MEGAsync/bio_comparadaII/result/Psychotria/final/output/Psychotria.new")
