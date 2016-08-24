# Creado por: Daniel Pabón
# email: dpabon@openmailbox.org

## ¿Cual prior hace sentido gamma o lognormal para el parametro shape de 
## la distribucion gama?

#######################################################################
##                          RECOMENDACIÓN                            ##
##        Utilizar una ventana de visualización externa              ##
##                              x11()                                ##
#######################################################################

## Primera parte del problema:

## Que distribución utilizar para recuperar el valor de shape, lognormal o gamma?

#comportamiento lognormal shape rate iguales
par(mfrow=c(2,1))
lo <-rlnorm(500, 0.1300, sdlog = 1)
plot(density(rgamma(n = 500, shape = lo[1], rate =lo[1])), ylim=c(0,4), xlim=c(0,4), main = "Recuperación a partir de distribución lognormal")
for(i in lo){
  lines(density(rgamma(n = 500, shape = i, rate = i)))
}
## muy pocos sitios evolucionan rapido, la mayoria de sitios evoluciona de forma lenta

#comportamiento gamma shape, rate iguales

lo <- rgamma(500, shape = 0.1300, rate = 0.1300)
plot(density(rgamma(n = 500, shape = lo[1], rate = lo[1]), from = 0, to = 4), ylim=c(0,4), xlim=c(0,80), main = "Recuperación a partir de distribución gamma")
for(i in lo){
  lines((density(rgamma(n = 500, shape = i, rate = i), from = 0, to = 4)))
}



## Al utilizar lognormal el prior es mas ajustado al valor de entrada
## mientras que al utilizar gamma  el espacio de muestreo es mas mayor
## al asignarle altas probabilidades a datos lejanos al valor inicial.
## Dado que el valor inicial fue recuperado con jmodeltest 
## hace mas sentido utilizar lognormal dado que al recuperar el valor de shape
## con gamma se terminan utilizando tasas de heterogeneidad 
## que no coinciden con la realidad.


## Segunda parte del problema:

## ¿Shape y rate diferentes, iguales o relacionados?


par(mfrow=c(2,1))
## lognormal con con rate = 2*shape
lo <-rlnorm(500, 0.1300, sdlog = 1)
plot(density(rgamma(n = 500, shape = lo[1], rate = lo[1]*2)), ylim=c(0,10), xlim=c(0,4), main = "Recuperación a partir de distribución lognormal")
for(i in lo){
  lines(density(rgamma(n = 500, shape = i, rate = i*2)))
}

## lognormal con con rate = 10*shape
lo <-rlnorm(500, 0.1300, sdlog = 1)
plot(density(rgamma(n = 500, shape = lo[1], rate = lo[1]*10), from = -2, to = 4), ylim=c(0,10), xlim=c(0,4), main = "Recuperación a partir de distribución lognormal")
for(i in lo){
  lines(density(rgamma(n = 500, shape = i, rate = i*10), from = -2, to = 4))
}


## lognormal con con rate = 1/shape
lo <-rlnorm(500, 0.1300, sdlog = 1)
plot(density(rgamma(n = 500, shape = lo[1], rate = 1/lo[1])), ylim=c(0,10), xlim=c(0,4), main = "Recuperación a partir de distribución lognormal")
for(i in lo){
  lines(density(rgamma(n = 500, shape = i, rate = 1/i), from = -2, to = 4))
}


## cuando el parametro rate es alto la escala disminuye inversamente por lo cual
## los valores se concentran mas entorno a la media (parametro de entrada).


## En conclusion:

## En mi caso donde poseo secuencias ITS donde pocos sitios cambian mucho
## hace mas sentido utilizar el shape reportado por jmodeltest con una desviación estandar
## de 1 para generar números aleatorios utilizando una distribución lognormal
## y luego pasar el valor como parametro para shape y rate = 1/shape
