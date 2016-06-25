library(ape)
library(seqinr)

## Números de acceso
datos <-read.csv("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos.csv")
directorios <- c("ITS_1", "ITS_2", "ITS_1+2")
## Descargando secuencias de Arctostaphylos

for(dic in 1:length(directorios)){
  setwd(paste("~/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/", directorios[dic], sep = ""))
  for(i in 1:length(datos[,dic+1])){
    if(is.na(datos[i,dic+1])==F){
      dato <- read.GenBank(datos[i,dic+1], species.names = T)
      write.dna(dato, paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), format = "fasta")
      #dato1 <- read.fasta(file =paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), seqtype = "DNA" )
      #write.fasta(dato1, names = paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""), file.out = paste(attr(dato, "specie"), "_", i, ".fasta", sep = ""))
    }
  }
}

prueba <- read.GenBank(datos[1,4])
prueba
class(datos[29,2])
is.na(datos[29,2])
