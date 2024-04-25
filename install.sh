#!/bin/bash
sudo apt update && sudo apt install -y git python3-pip python3-numpy python3-gdal python3-grib python3-netcdf4 ffmpeg
pip3 install matplotlib==2.1.1 scipy scikit_image folium moviepy similaritymeasures plotly Pillow scikit_learn scikit-image --user
#cupy==7.1.1 optionals

#clone the project
git clone http://pluton.hopto.org:8080/repositories/aedes_aegypti.git

#**Installing pybind11**
git clone https://github.com/pybind/pybind11.git
sudo cp pybind11/include/pybind11 /usr/include/  -r

#**Installing Eigen**
git clone https://gitlab.com/libeigen/eigen.git
sudo cp eigen/Eigen /usr/include/  -r

#**c++ binding**
cd aedes_aegypti && g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.8 src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so

#download weather data
wget  -r -l 1 -A .csv -np -nH -R index.html http://pluton.hopto.org:8080/data/public

#TODO: add line  for the cron

#add private folder
unzip ../private.zip -d data/

#to test
#python src/test 0;python src/test 1;python src/test 6;python src/test fit;python src/test plotFit;python src/test plotFitConf;
python3 src/tests.py 6
