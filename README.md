# Data

## IMERG

[LINK](https://search.earthdata.nasa.gov/)

GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree V07 (GPM_3IMERGDF) at GES DISC

- todo:
    * use the [API](https://gpm-api.readthedocs.io/en/latest/index.html) or [Google Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/NASA_GPM_L3_IMERG_V07#description)

## NCAR

[LINK](https://rda.ucar.edu/datasets/d083003/)

NCEP GDAS/FNL 0.25 Degree Global Tropospheric Analyses and Forecast Grids

- todo:
    * use the [API](https://github.com/NCAR/rda-apps-clients/tree/main)

These instructions assume that your are using [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or any of its flavours.

# Installing

## Clone the repository

```bash
git clone git@github.com:InstitutoGulich/aedes.git
```

## Create Conda environment 

```bash
conda env create --file environment.yml
```
## Create data folder and change its mode

```bash
mkdir -p data/public
chmod a+w data/public
```

Download the data tar file from [here](https://drive.google.com/file/d/1cUsIabnSyhezCHoGRMWhX3d-4qymQYrj/view?usp=sharing) and extract its content into the data/public folder.

```bash
tar xvf data-public_2024_09_24.tar.gz -C data/public/
```

## Compiling C++ library

```bash
g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I$CONDA_PREFIX/include/python3.8 -I$CONDA_PREFIX/include/eigen3 -I$CONDA_PREFIX/include/ src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so
```

In order to test this library, you can compile *main.cpp* in *src/cpp*:

```bash
g++ -std=c++17 -Wall -O3 -march=native -I$CONDA_PREFIX/include/python3.8 -I$CONDA_PREFIX/include -L./src/ src/cpp/main.cpp -o main.x
```

and then run 

```bash
./main.x
```
## Notebook

The notebook *example.ipynb* shows how you can use the library both *Python* and *C++* versions.


# Docker containers

```bash
docker-compose build
docker-compose up -d
```
