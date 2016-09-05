#!/bin/bash

## Script post-simulacion de arboles GEOSSE

cd ~/MEGAsync/bio_comparadaII/data/simulations/E1/
sed -i '2d' geo_*
echo "Escenario 1 LISTO!!"
cd ~/MEGAsync/bio_comparadaII/data/simulations/E2/
sed -i '2d' geo_*
echo "Escenario 2 LISTO!!"
