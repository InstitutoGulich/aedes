# Download data

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

# Conda environment 

```bash
conda env create --file environment.yml
```

# Compiling C++ library

```bash
g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I$CONDA_PREFIX/include/python3.8 -I$CONDA_PREFIX/include/ src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so
```

In order to test this library, you can compile main.cpp in src/cpp:

```bash
g++ -std=c++17 -Wall -O3 -march=native -I$CONDA_PREFIX/include/python3.8 -I$CONDA_PREFIX/include -L./src/ src/cpp/main.cpp -o main.x
```

and then run 

```bash
./main.x
```

# Dengue app

Create data folder

```bash
mkdir data
```

```bash
docker-compose build
docker-compose up -d
```
