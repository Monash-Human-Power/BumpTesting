# MHP_BumpTesting

This repo contains all code related to Bump Testing that was conducted from December 2017-March 2018.

- Acclerometers = contains code to run all four accelerometers which were used at different points in the testing.
- Filtering Code = contains matlab code to post-process the data to attempt to get some meaningful results.
- PIGPIO = contains software that needs to be installed to run the scripts on a new pi. To install:
cd pigpio
make -j4
sudo make install
