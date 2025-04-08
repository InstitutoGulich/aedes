import time
import numpy as np
import pandas as pd
import datetime

from utils import getPrecipitationsFromCsv
from config import Configuration
from model import Model

def myConf(dictionary=None):

    configuration = Configuration('example.cfg')

    if dictionary:
        for section in dictionary:
            for key in dictionary[section]:
                value = dictionary[section][key]
                if isinstance(value,list):
                    configuration.config_parser.set(section, key,','.join([str(x) for x in value]))
                else:
                    configuration.config_parser.set(section, key,str(value))

    configuration.validate()

    configuration.save('myConf.cfg')
    return configuration

if __name__ == '__main__':

    OUTPUT_FILENAME = 'output.csv'

    update = {'breeding_site':{ 
                                'amount':1
                            }
            }

    configuration = myConf(update)

    model = Model(configuration)
    t1 = time.time()
    time_range, results = model.solveEquations()
    t2 = time.time()
    print('Elapsed time: ', t2-t1)

    indexOf=lambda t: (np.abs(time_range-t)).argmin()

    start_datetime = datetime.datetime.strptime(configuration.getString('simulation','start_date'),'%Y-%m-%d')
    end_datetime = datetime.datetime.strptime(configuration.getString('simulation','end_date'),'%Y-%m-%d')
    dates = [(start_datetime + datetime.timedelta(days=t)) for t in time_range]

    parameters = model.parameters

    EGG    = parameters.EGG
    LARVAE = parameters.LARVAE
    PUPAE  = parameters.PUPAE
    ADULT1 = parameters.ADULT1
    ADULT2 = parameters.ADULT2
    WATER  = parameters.WATER
    OVIPOSITION = parameters.OVIPOSITION
    BS_a   = parameters.BS_a

    E = np.sum(results[:,EGG],axis=1)/BS_a
    L = np.sum(results[:,LARVAE],axis=1)/BS_a
    A = (results[:,ADULT1]+results[:,ADULT2])/BS_a

    lwO = np.array([results[indexOf(t),OVIPOSITION] - results[indexOf(t-7),OVIPOSITION] for t in time_range])
    lwO_mean = np.array([lwO[indexOf(t-7):indexOf(t+7)].mean(axis=0) for t in time_range])
    O = np.sum(lwO_mean,axis=1)/BS_a

    T  = parameters.weather.T
    RH = parameters.weather.RH
    P  = parameters.weather.p

    location = parameters.location['name']

    T = T(time_range) - 273.15
    P = P(time_range)
    RH = RH(time_range)

    df = pd.DataFrame({'date':dates,'E':E,'L':L,'A':A,'O':O,'p':P,'T':T,'RH':RH})
    df.set_index('date',inplace=True)

    df.to_csv(OUTPUT_FILENAME,index=True)