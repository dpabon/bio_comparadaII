## Escalando revBayes

#Tiempos de salida

#mpi + sequential =
#real	466m32.776s
#user	773m12.036s
#sys	0m9.056s

#mpi =
#real	165m29.557s
#user	659m28.720s
#sys	0m13.172s


#sequential = 
#real	620m4.992s
#user	619m56.872s
#sys	0m4.636s

tiempos <- c(620, 466, 165)
tiempos <- tiempos/60
names(tiempos) <-c("SEQ", "MPI+SEQ", "MPI" )
barplot(tiempos, ylab = "Horas", col = "#f44336", border = F, 
        main = "Tiempo de ejecución \n 10 millones de generaciones
        y arbol de mayor probabilidad")

