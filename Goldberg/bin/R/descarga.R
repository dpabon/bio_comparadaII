datos <-read.csv("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos.csv")
setwd("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/ITS_1/")
library(ape)

for(i in 1:length(datos[,2])){
  if(i==NULL){
    print("no hay secuencia disponible")
  }else{
    dato <- read.GenBank(i, species.names = T)
    write.dna(dato, paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), format = "fasta")
  }
}
