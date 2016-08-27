#!/bin/bash

## Script prioris flat Arctostaphylos

echo 'Instalando compiladores...'
sudo apt-get install build-essential
sudo apt-get install linux-headers-`uname -r`
sudo apt-get install checkinstall make automake cmake autoconf
echo 'Compiladores instalados'

## Programas de utilidad
sudo apt-get install git mpich unzip vim libboost-all-dev libboost-date-time1.58.0 libboost-filesystem1.58.0 libboost-iostreams1.58.0 libboost-python1.58.0 libboost-serialization1.58.0 libboost-system1.58.0 libboost-thread1.58.0 libboost-atomic-dev libboost-atomic1.58.0 libboost-chrono-dev libboost-chrono1.58.0 libboost-thread-dev libboost-dev libboost1.58-dev

sudo reboot
### Instalando revbayes 1.0.1
git clone https://github.com/revbayes/revbayes.git revbayes-master
cd revbayes-master/projects/cmake/
./build.sh -mpi true
sudo mv rb-mpi /usr/local/bin/
## Archivos de corrida
cd ~
mkdir MEGAsync/
cd MEGAsync/
git clone https://github.com/dpabon/bio_comparadaII
mkdir flat
cd flat

nohup mpirun -np 4 rb-mpi /home/dane/MEGAsync/bio_comparadaII/data/Arctostaphylos/revbayes/Arctostaphylos_clock.rev &

echo 'Arctostaphylos flat corriendo...'
echo 'Vuelva mas tarde'
