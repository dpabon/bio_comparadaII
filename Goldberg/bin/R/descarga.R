library(ape)
## Numeros de acceso
datos <-read.csv("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos.csv")
directorios <- c("ITS_1", "ITS_2", "ITS_1+2")

## Descargando las secuencias
for(dic in 1:length(directorios)){
  setwd(paste("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/", directorios[dic], sep = ""))
  for(i in 1:length(datos[,dic+1])){
    if(is.null(i)==T){
      print("no hay secuencia disponible")
    }else{
      dato <- read.GenBank(datos[i,dic+1], species.names = T)
      write.dna(dato, paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), format = "fasta")
    }
  }
}


for(i in 1:length(datos[,2])){
  if(is.null(i)==T){
    print("no hay secuencia disponible")
  }else{
    dato <- read.GenBank(datos[i,2], species.names = T)
    write.dna(dato, paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), format = "fasta")
  }
}

setwd("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/ITS_2/")

for(i in 1:length(datos[,3])){
  if(is.null(i)==T){
    print("no hay secuencia disponible")
  }else{
    dato <- read.GenBank(datos[i,3], species.names = T)
    write.dna(dato, paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), format = "fasta")
  }
}

