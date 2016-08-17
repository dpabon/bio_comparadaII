#!/bin/bash

## dependencias necesarias
echo 'Instalando compiladores...'
sudo apt-get install build-essential
sudo apt-get install linux-headers-`uname -r`
sudo apt-get install checkinstall make automake cmake autoconf
echo 'Compiladores instalados'

## Programas de utilidad
sudo apt-get install git unzip vim

### Instalando revbayes 1.0.0
wget https://github.com/revbayes/revbayes/archive/v1.0.0-release.zip
unzip v1.0.0-release.zip
cd revbayes-1.0.0-release/projects/cmake/
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
