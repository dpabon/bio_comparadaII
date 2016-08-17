#!/bin/bash

cd /home/user/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/ITS_1
rm compendio_ITS_1.fas
cat * > compendio_ITS_1.fas
cd /home/user/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/ITS_2
rm compendio_ITS_2.fas
cat * > compendio_ITS_2.fas
cd /home/user/MEGAsync/bio_comparadaII/Goldberg/data/Arctostaphylos/ITS_1+2
rm compendio_ITS_1+2.fas
cat * > compendio_ITS_1+2.fas
