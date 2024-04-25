import os
import tempfile
import datetime
import numpy as np
from config import Configuration

try:
    from otero_precipitation_wrapper import Model as _Model
except ImportError:
    os.system('g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.8 src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so')
    from otero_precipitation_wrapper import ModelWrapper as _Model

class ParametersDecorator:
    def __init__(self,parameters):
        self.parameters=parameters
    def __getattr__(self, name):
        attribute=getattr(self.parameters,name)
        if(name in ['EGG','LARVAE','PUPAE','WATER','OVIPOSITION']):
            return slice(attribute.first(),attribute.first()+attribute.size())
        elif(name in ['location']):
            return {'name':attribute}
        else:
            return attribute

class Model:
    def __init__(self, configuration=Configuration('resources/otero_precipitation.cfg')):
        self.config_filename=tempfile.NamedTemporaryFile(suffix='.cfg').name
        with open(self.config_filename, 'w') as configfile:
            configuration.config_parser.write(configfile)

        self._model=_Model(self.config_filename)
        self.start_date=datetime.datetime.strptime(self._model.start_date,'%Y-%m-%d').date()
        self.end_date=datetime.datetime.strptime(self._model.end_date,'%Y-%m-%d').date()
        self.time_range=np.array(self._model.time_range)
        self.parameters=ParametersDecorator(self._model.parameters)
        self.warnings=['warnings not implemented']


    def solveEquations(self):
        self._model.solveEquations()
        self.Y=np.array(self._model.Y)
        return self.time_range,self.Y
