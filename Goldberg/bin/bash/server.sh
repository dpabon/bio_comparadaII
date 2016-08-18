#!/bin/bash

## dependencias necesarias
echo 'Instalando compiladores...'
sudo apt-get install build-essential
sudo apt-get install linux-headers-`uname -r`
sudo apt-get install checkinstall make automake cmake autoconf
echo 'Compiladores instalados'

## Programas de utilidad
sudo apt-get install git unzip vim libboost-date-time1.58.0 libboost-filesystem1.58.0 libboost-iostreams1.58.0 libboost-python1.58.0 libboost-serialization1.58.0 libboost-system1.58.0 libboost-thread1.58.0

### Instalando revbayes 1.0.0
git clone https://github.com/revbayes/revbayes.git revbayes-master
cd revbayes-master/projects/cmake/
./build.sh
sudo mv rb /usr/local/bin/
## Archivos de corrida

mkdir MEGAsync/
cd MEGAsync/
git clone https://github.com/dpabon/bio_comparadaII
mkdir results/
cd results/

nohup rb /home/dane/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/revbayes/Arctostaphylos_clock.rev &

nohup rb /home/dane/MEGAsync/bio_comparadaII/Goldber/data/Ceanothus/revbayes/Ceanothus_clock.rev &

echo "Arctostaphylos y Ceanothus corriendo vuelva mas tarde"
