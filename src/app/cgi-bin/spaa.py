import os
import json
import datetime
import numpy as np
import sys
os.chdir('../../../')
sys.path.append('./src')
import utils
from config import Configuration
from otero_precipitation_wrapper_wrapper import Model
from equations import diff_eqs

SCENARIO_FULLNAMES={'rc':'condiciones_regulares','it':'temperatura_aumentada','dt':'temperatura_disminuida','ip':'precipitacion_aumentada','dp':'precipitacion_disminuida'}

def daysSinceEpoch(start_datetime,days):
    epoch = datetime.datetime.utcfromtimestamp(0)
    return (start_datetime + datetime.timedelta(days=days) - epoch).total_seconds() * 1000.0

def applyScenario(location,scenario):
    if(scenario in ['it','dt']):
        if(scenario=='it'):#increased temperature
            a=1.2
        elif(scenario=='dt'):#decreased temperature
            a=0.8
        location=utils.scaleWeatherConditions('data/public/',location,2,a)

    if(scenario in ['ip','dp']):
        if(scenario=='ip'):#increased precipitation
            a=1.2
        elif(scenario=='dp'):#decreased precipitation
            a=0.8
        location=utils.scaleWeatherConditions('data/public/',location,4,a)

    return location

def runSimulation(GET):
    start_datetime=datetime.datetime.strptime(GET.get('start_date'),'%Y-%m-%d')
    start_date=start_datetime.date()
    end_date=datetime.datetime.strptime(GET.get('end_date'),'%Y-%m-%d').date()
    location=GET.get('location')
    scenario=GET.get('scenario')
    location=applyScenario(location+'.full', scenario)
    configuration=Configuration('resources/1c.cfg',
        {'simulation':{
            'start_date':start_date,
            'end_date':end_date,
            },
        'location':{
            'name':str(location)#unicode to str
        },
        })
    model=Model(configuration)

    time_range,Y= model.solveEquations()
    EGG,LARVAE,ADULT1,ADULT2,OVIPOSITION=model.parameters.EGG,model.parameters.LARVAE,model.parameters.ADULT1,model.parameters.ADULT2,model.parameters.OVIPOSITION
    BS_a=configuration.getFloat('breeding_site','amount')
    indexOf=lambda t: (np.abs(time_range-t)).argmin()

    E=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,EGG])/BS_a ] for i,t in enumerate(time_range)]
    L=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,LARVAE])/BS_a ] for i,t in enumerate(time_range)]
    A=[ [ daysSinceEpoch(start_datetime,t), (Y[i,ADULT1]+Y[i,ADULT2])/BS_a ] for i,t in enumerate(time_range)]

    lwO=np.array([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range])
    lwO_mean=np.array([lwO[indexOf(t-7):indexOf(t+7)].mean(axis=0) for t in time_range])
    O=[ [ daysSinceEpoch(start_datetime,t), np.sum(lwO_mean[i])/BS_a ] for i,t in enumerate(time_range)]

    weather=model.parameters.weather
    precipitations = utils.getPrecipitationsFromCsv('data/public/'+location+'.csv',start_date,end_date)
    p=[ [ daysSinceEpoch(start_datetime,t), precipitations[int(t)] ] for t in time_range]
    T=[ [ daysSinceEpoch(start_datetime,t), float(weather.T(t)) - 273.15 ] for t in time_range]
    RH=[ [ daysSinceEpoch(start_datetime,t), float(weather.RH(t)) ] for t in time_range]

    return json.dumps({
                        'population':[{'name':'Huevos','data':E,'type':'scatter'},{'name':'Oviposicion','data':O,'type':'scatter'}],
                        'population_II':[{'name':'Larvas','data':L,'type':'scatter'},{'name':'Adultos','data':A,'type':'scatter'}],
                        'weather':[{'name':'Temperatura','data':T,'type':'scatter'},{'name':'Humedad Relativa','data':RH,'type':'scatter'},{'name':'Precipitacion','data':p,'type':'bar'}]
                        })
