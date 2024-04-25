#coding: utf-8
import os
import sys
import utils
import datetime
import numpy as np
import multiprocessing as mp
from config import Configuration
from configparser import ConfigParser
from otero_precipitation import Model
from equations import diff_eqs,vR_D
import equations
from spatial_equations import diff_eqs as spatial_diff_eqs
import similaritymeasures as sm
import plotly.graph_objs as go
from plotly import tools
import equation_fitter
import unidecode

def runSpatial():
    configuration=Configuration('resources/otero_precipitation.cfg')
    configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
    model=Model(configuration)
    #modify some parameters to make them spatial
    parameters=model.parameters
    n=parameters.n
    parameters.FLYER=3*n+2#in R
    parameters.diff=configuration.getFloat('biology','diffusion')#diffusion-like coefficient
    parameters.P,warning=utils.getPreferenceMatrix('data/public/goodness/'+configuration.getString('biology','goodness'),patch_size=100)
    model.warnings.append(warning)
    HEIGHT,WIDTH=parameters.P.shape[:2]
    parameters.initial_condition=np.append(parameters.initial_condition,[0])#append flyers
    parameters.initial_condition=(parameters.initial_condition*utils.getY0FactorMatrix(HEIGHT,WIDTH)[:,:,np.newaxis]).reshape(HEIGHT*WIDTH*(3*n+3))#TODO:find a better way of introducing initial conditions to spatial
    parameters.vBS_a=parameters.BS_a*np.ones((HEIGHT,WIDTH))#TODO:Estimate BS_a as in section 8.3 Otero 2008, "Estimation of the breeding site density"
    parameters.vBS_d=parameters.vBS_d*np.ones((HEIGHT,WIDTH,n))#ASSUMPTION: distribuition is constant on (x,y)
    parameters.vAlpha=(parameters.vAlpha0*np.ones((HEIGHT,WIDTH,n)) )/parameters.vBS_a[:,:,np.newaxis]
    theta=parameters.vBS_a/150.
    tdep=0.229#Average time for egg deposition Christophers(1960)
    parameters.ovr= np.where(parameters.vBS_a<=150, theta/tdep, 1/tdep)

    for warning in model.warnings:
        print('# WARNING: ' + warning)

    #solve the equations
    time_range,Y=model.solveEquations(equations=utils.ProgressEquations(model,spatial_diff_eqs),method='cuda_rk' )
    Y=Y.reshape(Y.shape[0],HEIGHT,WIDTH,3*n + 3)
    np.save('out/Y.npy',Y)
    #time_range,Y=model.time_range,np.load('out/Y.npy')#to debug video

    EGG,LARVAE,PUPAE,ADULT1,FLYER,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.FLYER,parameters.ADULT2
    stages={'E':EGG, 'A':[ADULT1,FLYER,ADULT2]}
    for key in stages:
        print('Creating animation for %s...'%key)
        matrix=np.sum(Y[:,:,:,stages[key]],axis=3)
        matrix=matrix/matrix.max()
        start_date=configuration.getDate('simulation','start_date')
        getTitle=lambda i: datetime.timedelta(days=time_range[i])+start_date
        utils.createAnimation('out/%s'%key,matrix,getTitle,time_range.max())# 1 day : 1 second

from scipy import stats
def calculateMetrics(time_range,lwO_mean,ovitrap_eggs_i):
    valid_ovi_idx=~np.isnan(ovitrap_eggs_i)
    rmse=utils.rmse(ovitrap_eggs_i[valid_ovi_idx],lwO_mean[valid_ovi_idx])
    cort=utils.cort(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])
    pearson,p_value=stats.pearsonr(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])

    reversed_valid_ovi_idx=valid_ovi_idx[::-1]
    first,last=np.argmax(valid_ovi_idx), len(reversed_valid_ovi_idx)-np.argmax(reversed_valid_ovi_idx)-1
    x=np.array([[time_range[idx],lwO_mean[idx]] for idx in range(first,last)])
    y=np.array([ [time_range[idx],ovitrap_eggs_i[idx] ] for idx,isValid in enumerate(valid_ovi_idx) if isValid])
    fd=sm.frechet_dist(x,y)
    dtw, path = sm.dtw(x, y)
    D_1=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=1)
    D=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])
    D_2=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=2)
    D_4=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=4)

    return rmse, cort,pearson,fd,dtw,D,D_1,D_2,D_4

def runCases(case):
    if(case==0):
        ovi_range=range(1,151)
        errors_by_height=np.empty((15,151,9))
        errors_by_height[:]=np.nan
        for h in range(1,15):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')#TODO:fix data and
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(20)))# uncomment these two
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,Y=model.solveEquations()

            #errors=[[1e15,-1,-1,1e15,1e15,1e15]]*151#just to fill the ovitrap 0 that do not exist in reality
            for ovitrap_id in ovi_range:
                values=utils.getOvitrapEggsFromCsv('data/private/ovitrampas_2017-2019.full.csv' ,ovitrap_id)
                ovi=np.array(equation_fitter.populate(model.time_range,model.start_date,values))
                ovi=np.array(ovi,dtype=np.float)#this change None for np.nan

                indexOf=lambda t: (np.abs(time_range-t)).argmin()
                OVIPOSITION=model.parameters.OVIPOSITION
                BS_a=model.parameters.BS_a
                O=Y[:,OVIPOSITION]
                lwO=np.sum([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range],axis=1)/BS_a#calculate the difference,sum up all levels, and divide by amount of containers

                errors_by_height[h][ovitrap_id]=calculateMetrics(time_range,lwO,ovi)


            errors=errors_by_height[h]
            rmse, cort,pearson,fd,dtw,D,D_1,D_2,D_4=[errors[:,i] for i in range(errors.shape[1])]
            print('''h: %scm.
                         id,    score
                rmse:    %3s,   %s
                cort:    %3s,   %s
                pearson: %3s,   %s
                fd:      %3s,   %s
                dtw:     %3s,   %s
                D:       %3s,   %s'''%
                (h,
                np.nanargmin(rmse),np.nanmin(rmse),
                np.nanargmin(cort),np.nanmin(cort),
                np.nanargmin(pearson),np.nanmin(pearson),
                np.nanargmin(fd),np.nanmin(fd),
                np.nanargmin(dtw),np.nanmin(dtw),
                np.nanargmin(D),np.nanmin(D)
                ) )

            #print first N
            for i in range(errors.shape[1]):
                f=errors[:,i]
                ovi_sorted=np.argsort(f)
                N=5
                print('''h: %scm.
                     id,    score
                     f:    %3s  <---->  %s'''%
                     (h,
                     ', '.join(map(str,ovi_sorted[:N])), ', '.join(map(str,f[ovi_sorted[:N]]))
                     ))
                #print last N(excluding the ficticius one)
                print('''h: %scm.
                     id,    score
                     f:    %3s  <---->  %s'''%
                     (h,
                     ', '.join(map(str,ovi_sorted[-N-1:-1])), ', '.join(map(str,f[ovi_sorted[-N-1:-1]]))
                     ))

        np.save('out/errors_by_height.npy',errors_by_height)
        print(errors_by_height.shape)

    if(case==1):
        h=10.
        configuration=Configuration('resources/1c.cfg')
        configuration.config_parser.set('location','name','cordoba.full')
        location=configuration.getString('location','name')
        configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
        configuration.config_parser.set('breeding_site','height',str(h))
        model=Model(configuration)
        time_range,Y=model.solveEquations()
        utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list([151]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
        title=location.replace('.full','').replace('_',' ').title(),
        xaxis_title='Date',
        yaxis_title='Eggs')

        #utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.(Oct-Nov-Dic just prom available)'%(mf,h))
        #utils.showPlot(utils.plot(model,subplots=[{'pa':''}]))
        print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
        print(model.warnings)

        #is OEquations perturbing the result somehow?No, the results match.
        #model2=Model(configuration)
        #time_range2,initial_condition2,Y2=model2.solveEquations(method='rk')
        #print(np.linalg.norm((Y[:,:model.parameters.OVIPOSITION.start]-Y2)))

    if(case==2):
        for mf  in [0.,3.]:
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','manually_filled',str(mf))
            h=configuration.getArray('breeding_site','height')[0]
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.'%(mf,h))
            print(model.warnings)

    if(case==3):
        for year in range(2016,2019):
            for month in range(1,13):
                start_date,end_date=datetime.date(year,month,1),datetime.date(year+int(month/12),month%12 +1,1)
                precipitations = utils.getPrecipitationsFromCsv(sys.argv[2],start_date,end_date)
                print('period: %s to %s        %smm.'%(start_date,end_date,np.sum(precipitations)))
    if(case==4):
        ovi_range=range(1,151)
        ovi_mean=[1e10]*151
        for ovitrap_id in ovi_range:
            OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.full.csv'
            values=utils.getOvitrapEggsFromCsv(OVITRAP_FILENAME,ovitrap_id)
            dates=values.keys()
            ovi_a=[values[date][0] if date in values else None for date in dates]#TODO:WARNING!this will repeat values if model granularity is not 1 value per day.
            ovi_a=np.array(ovi_a,dtype=np.float)#this change None for np.nan
            ovi_b=[values[date][1] if (date in values and len(values[date])>1) else None for date in dates]
            ovi_b=np.array(ovi_b,dtype=np.float)#this change None for np.nan
            ovi_mean[ovitrap_id]=np.nanmean(np.abs(ovi_a-ovi_b)/(ovi_a+ovi_b))

        ovi_mean=np.array(ovi_mean)
        ovi_ordered=np.argsort(ovi_mean)
        print(ovi_ordered)
        print(ovi_mean[ovi_ordered])
    if(case==5):
        for location in ['cordoba.full.weather-2019-01-01','cordoba.full.weather-2019-02-01']:
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location)
            start_date,end_date=utils.getStartEndDates('data/public/'+location+'.csv')
            configuration.config_parser.set('simulation','end_date',str(end_date))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list([143]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title=location.replace('.full.weather-',' ').title(),
            xaxis_title='Date',
            yaxis_title='Number of eggs')
            print(model.warnings)

    if(case==6):
        config_parser = ConfigParser()
        config_parser.read('resources/get_weather.cfg')
        for location in ['cordoba']:#config_parser.sections():
            h=10.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location+'.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            utils.showPlot(utils.plot(model,subplots=[{'cd':'','A1+A2':'','W':'','f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title=location.replace('_',' ').title(),
            xaxis_title='Date',
            yaxis_title='')
            print(model.warnings)
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))

    if(case==7):
        configuration=Configuration('resources/1c.cfg')
        configuration.config_parser.set('location','name','cordoba.full')
        configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
        model=Model(configuration)
        time_range,Y=model.solveEquations()
        utils.showPlot(utils.plot(model,subplots=[{'E':':','A1+A2':'','f':[utils.safeAdd]}]),
            title='Cordoba',
            xaxis_title='Date',
            yaxis_title='Individuals')
        print(model.warnings)

    if(case==8):
        LOCATION='cordoba'
        PLOT_START_DATE=datetime.date(2018,12,15)#maybe by parameter?
        FORECAST=51
        PLOT_END_DATE=PLOT_START_DATE+datetime.timedelta(FORECAST)
        DATA_FOLDER='data/public/'
        HISTORY_FOLDER=DATA_FOLDER  + '.history/'
        data_A,data_O,data_W=[],[],[]
        filenames=os.listdir(HISTORY_FOLDER)
        filenames.sort()
        i=-1
        for filename in  filenames:
            if(filename=='.empty' or not filename.startswith(LOCATION)): continue
            location,year,month,day=filename.replace('.full.weather','').replace('.csv','').split('-')
            simulation_date=datetime.date(int(year),int(month),int(day))
            if( not (PLOT_START_DATE<=simulation_date<=PLOT_END_DATE)): continue#the equals is because we want to  have one last curve with no forecast
            i=i+1
            if(i==0):
                color='rgb(0, 0, 0)'
            elif(i==FORECAST-1):
                color='rgb(255, 0, 255)'
            else:
                color = 'rgb(%s, %s, 0)'%(int( (1-i/FORECAST) * 255),int(i/FORECAST * 255))
            if(i%7!=0 and i!=FORECAST-1): continue#every 7 days but the last one must be shown.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','.history/'+filename.replace('.csv',''))
            configuration.config_parser.set('simulation','end_date',str(PLOT_END_DATE))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            data_A+=utils.plot(model,subplots=[{'A1+A2':str(simulation_date),'f':[utils.safeAdd]}],plot_start_date=PLOT_START_DATE,color=color)
            data_O+=utils.plot(model,subplots=[{'lwO':str(simulation_date)  ,'f':[utils.safeAdd]}],plot_start_date=PLOT_START_DATE,color=color)
            data_W+=utils.plot(model,subplots=[{'W':str(simulation_date)    ,'f':[]             }],plot_start_date=PLOT_START_DATE,color=color)
            print(filename,model.warnings)


        x=[]
        y=[]
        for serie in data_A:
            x+= [( datetime.datetime.strptime(data_A[-1]['name'],'%Y-%m-%d') - datetime.datetime.strptime(serie['name'],'%Y-%m-%d') ).days]#TODO:we depend on using simulation date as name
            S_last,S_i=np.array(data_A[-1]['y']),np.array(serie['y'])
            y+= [np.sum(np.abs(S_last-S_i)/S_last)*1/len(S_i) ]#relative difference mean. 0<=||S_last-S_i||/(S_last+S_i) <= 1 by triangular inequality. => sum(...)/n in [0,1].
            assert len(S_i)==len(S_last)
        x,y=[e for e in reversed(x)],[e for e in reversed(y)]

        utils.showPlot(data_A,title='Adults in '+LOCATION.title(),xaxis_title='Date',yaxis_title='Individuals')
        utils.showPlot(data_O,title='Oviposition in '+LOCATION.title(),xaxis_title='Date',yaxis_title='Eggs')
        utils.showPlot(data_W,title='Water in '+LOCATION.title(),xaxis_title='Date',yaxis_title='cm.')
        utils.showPlot([go.Scatter(x=x,y=y, name='')],title='',xaxis_title='Days',yaxis_title=u'\u03B4')

    names=['rmse', 'cort','pearson','fd','dtw','D','D_1','D_2','D_4']
    d=5
    if(case==9):
        errors_by_height=np.load('out/errors_by_height.npy')
        for d,name in enumerate(names):
            utils.showPlot([go.Surface(z=errors_by_height[:,:,d] )],title=name,scene=dict(xaxis=dict(title='Ovitrap id'),yaxis=dict(title='Height')) )
    if(case==10):
        errors_by_height=np.load('out/errors_by_height.npy')
        for h in range(1,15):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            o_h=[]
            for i in range(1,151):
                if np.nanargmin(errors_by_height[:,i,d])==h: o_h+=[i]

            if(o_h):
                utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list(o_h),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
                    title=names[d]+' Height: %scm.'%h,
                    xaxis_title='Fecha',
                    yaxis_title='Nº de huevos')

    if(case==11):
        errors_by_height=np.load('out/errors_by_height.npy')
        for d,name in enumerate(['D']):
            heights=[1,3,6,8]
            fig = tools.make_subplots(rows=len(heights), cols=1)
            for row,h in enumerate(heights):#range(1,errors_by_height.shape[0]):
                fig.append_trace(go.Scatter(x=np.array(range(0,errors_by_height.shape[1])), y=errors_by_height[h,:,d],name='%scm.'%h, mode='markers'),row+1,1)
                #fig['layout']['yaxis'+str(h)].update(range=[1,np.nanmax(errors_by_height[:,:,d])])
            fig['layout']['title'] = ''#https://community.plot.ly/t/subplots-how-to-add-master-axis-titles/13927
            fig['layout']['annotations']=[go.layout.Annotation(x=0.5,y=-0.15,showarrow=False,text="Ovitrap identifier",xref="paper",yref="paper",),
                                            go.layout.Annotation(x=-0.07,y=0.5,showarrow=False,text="Dissimilarity Index  D",textangle=-90,xref="paper",yref="paper")]
            fig['layout']['font'] = dict(family="Courier New, monospace",size=24,color="#090909")
            fig['layout']['autosize']=True
            fig['layout']['margin']=dict(l=120,b=120)

            utils.showPlot(fig)

    if(case==12):
        heights=[1,3,6,8,10]
        ovis=[13,134,54,122,151]
        for i,h in enumerate(heights):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            utils.showPlot(utils.plot(model,subplots=[{'lwO':'','O':list([ovis[i]]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title='%scm.'%h,
            xaxis_title='Date',
            yaxis_title='Eggs')
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
            print(model.warnings)


    if(case==13):
        errors_by_height=np.load('out/errors_by_height.npy')
        a=[-1]+[np.nanmin(errors_by_height[:,o,5]) for o in range(1,151)]
        idx=np.argsort(a)
        print(idx[-5:])

    if(case==14):
        config_parser = ConfigParser()
        config_parser.read('resources/get_weather.cfg')
        for location in ['bahia_blanca','general_roca','cordoba','tartagal','santa_fe']:#config_parser.sections():
            h=10.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location+'.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,Y=model.solveEquations()

            utils.showPlot(utils.plot(model,subplots=[{'lwO':'','f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title=location.replace('_',' ').title(),
            xaxis_title='Date',
            yaxis_title='Eggs')
            print(model.warnings)
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))

    if(case==15):
        errors_by_height=np.load('out/errors_by_height.npy')
        matrix=errors_by_height[:,:,3:6]
        matrix/=np.nanmax(matrix,axis=(0,1))#normalize
        plt.imshow(np.uint8(matrix*255),interpolation='nearest', aspect='auto')
        plt.show()
    if(case==16):
        heights=[2,5,10,15]
        fig = tools.make_subplots(rows=2, cols=2,subplot_titles=['%scm.'%h for h in heights])
        for i,h in enumerate(heights):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            model.solveEquations()
            traces=utils.plot(model,subplots=[{'PO':''}],plot_start_date=datetime.date(2017,10,1))
            for trace in traces:
                fig.append_trace(trace,int(i/2) +1,i%2 +1)
        utils.showPlot(fig)
    if(case==17):
        errors_by_height=np.load('out/errors_by_height.npy')
        d=5
        utils.showPlot([go.Heatmap(z=errors_by_height[:,:,d])])

    if(case==18):
        config_parser = ConfigParser()
        config_parser.read('resources/get_weather.cfg')
        for stage,title in {'E':'Eggs population','L':'Larvae population','P':'Pupae population','A1+A2':'Adults population'}.items():
            h=10.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','height',str(h))
            configuration.config_parser.set('breeding_site','amount','1')
            model=Model(configuration)
            time_range,Y=model.solveEquations()

            utils.showPlot(utils.plot(model,subplots=[{stage:'','f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title='',
            xaxis_title='Date',
            yaxis_title=title)
            print(model.warnings)
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
    if(case==19):
        MAX=50000
        model=Y=time_range=None
        for i in range(1,MAX):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            print('\r %s%%'%round(i/MAX * 100,2), end='')

    if(case==20):
        dataO,dataA=[],[]
        for location in ['balnearia','marull','rio_cuarto','unquillo','san_francisco','cordoba']:
            h=10.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location+'.full')
            location=configuration.getString('location','name')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
            configuration.config_parser.set('breeding_site','height',str(h))
            configuration.config_parser.set('breeding_site','amount','1')
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            location_label=location.replace('.full','').replace('_',' ').title()
            dataO+=utils.plot(model,subplots=[{'lwO':location_label,'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))
            dataA+=utils.plot(model,subplots=[{'A1+A2':location_label,'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
            print(model.warnings)
        dataO+=utils.plot(model,subplots=[{'cd':''}])
        dataA+=utils.plot(model,subplots=[{'cd':''}])
        utils.showPlot(dataO,title='Oviposición',xaxis_title='Fecha',yaxis_title='Huevos')
        utils.showPlot(dataA,title='Hembra Adulta por criadero',xaxis_title='Fecha',yaxis_title='Adultos')

    if(case==21):
        dataO,dataA=[],[]
        lines=[line.strip().split(',') for line in open('data/private/Registros de presencia-ausencia de Ae. aegypti (por bibliografía y muestreados).csv').readlines()[1:]]
        for line in lines:
            h=10.
            location='%s'%unidecode.unidecode(line[3]).lower().strip().replace(' ','_')
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location)#+'.full'
            location=configuration.getString('location','name')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(-1)))
            configuration.config_parser.set('breeding_site','height',str(h))
            configuration.config_parser.set('breeding_site','amount','1')
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            location_label=location.replace('.full','').replace('_',' ').title() + ' (%s)'%line[-2] + ('*' if line[-1] else '')
            dataO+=utils.plot(model,subplots=[{'lwO':location_label,'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))
            dataA+=utils.plot(model,subplots=[{'A1+A2':location_label,'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))
            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
            print(model.warnings)
        dataO+=utils.plot(model,subplots=[{'cd':''}])
        dataA+=utils.plot(model,subplots=[{'cd':''}])
        utils.showPlot(dataO,title='Oviposición',xaxis_title='Fecha',yaxis_title='Huevos')
        utils.showPlot(dataA,title='Hembra Adulta por criadero',xaxis_title='Fecha',yaxis_title='Adultos')

    if(case==22):#generalized 21
        data={'lwO':[],'L':[],'P':[],'A1+A2':[],'W':[],'T':[],'pa':[],'RH':[],'clc':[],'LEx':[],'sP':[],'Lr0':[]}
        lines=[line.strip().split(',') for line in open('data/private/Registros de presencia-ausencia de Ae. aegypti (por bibliografía y muestreados).csv').readlines()[1:]]
        for line in lines:
            if(not line[-2]=='+' or line[-1]): continue
            h=10.
            location='%s'%unidecode.unidecode(line[3]).lower().strip().replace(' ','_')
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location)#+'.full'
            location=configuration.getString('location','name')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(-1)))
            configuration.config_parser.set('breeding_site','height',str(h))
            configuration.config_parser.set('breeding_site','amount','1')
            model=Model(configuration)
            time_range,Y=model.solveEquations()
            location_label=location.replace('.full','').replace('_',' ').title() + ' (%s)'%line[-2] + ('*' if line[-1] else '')
            for key in data:
                data[key]+=utils.plot(model,subplots=[{key:location_label,'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))

            print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
            print(model.warnings)

        #add mean curves
        for key in data:
            y=[]
            for serie in data[key]: y+=[serie['y']]
            data[key]+=[go.Scatter(x=data[key][0]['x'],y=np.mean(y,axis=0), name=key+' mean' )]
            data[key]+=[go.Scatter(x=data[key][0]['x'],y=np.mean(y,axis=0)+np.std(y,axis=0), name=key+' mean + std' )]
            data[key]+=[go.Scatter(x=data[key][0]['x'],y=np.mean(y,axis=0)-np.std(y,axis=0), name=key+' mean - std' )]

        exportCSV(data)

        for key in data:
           data[key]+=utils.plot(model,subplots=[{'cd':''}])
           utils.showPlot(data[key],title=key,xaxis_title='Fecha',yaxis_title='')

    if(case==23):
        lwO_values=[]
        A_values=[]
        O_PLUS_STD=260.
        A_MAX=35./15.
        filenames=os.listdir('data/public/XY/')
        for i,filename in enumerate(filenames):
            if(not filename.startswith('XY') or 'full' not in filename): continue
            h=10.
            BS_a=1
            #>sed -i 's/--/0.015/g' *.csv#sanitize files
            print('%s (%s/%s)'%(filename,i,len(filenames)),end='\r')#
            location='XY/'+filename.replace('.csv','')
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location)#+'.full'
            location=configuration.getString('location','name')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(-1)))
            configuration.config_parser.set('breeding_site','height',str(h))
            configuration.config_parser.set('breeding_site','amount',str(BS_a))
            model=Model(configuration)
            BS_a,vBS_d,m,n,OVIPOSITION,ADULT1,ADULT2=model.parameters.BS_a,model.parameters.vBS_d,model.parameters.m,model.parameters.n,model.parameters.OVIPOSITION,model.parameters.ADULT1,model.parameters.ADULT2
            time_range,Y=model.solveEquations()
            indexOf=lambda t: (np.abs(time_range-t)).argmin()
            lwO=np.array([ (Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION]).reshape(m,n).sum(axis=0) for t in time_range])/(BS_a*vBS_d)#here index and day number are equal
            lwO=lwO[:,0]#if multiple container w assume the first one is ovitrap and the rest are wild containers
            A=(Y[:,ADULT1]+Y[:,ADULT2])/BS_a

            trash,lat,lon=filename.replace('.full.csv','').split('_')
            lwO_values.append([[lat,lon,lwO_d/O_PLUS_STD + 1e-4] for lwO_d in lwO ])#d in days#for some reason weght 0 is not allowed and produce an spasmodic result
            A_values.append([[lat,lon,A_d/A_MAX + 1e-4] for A_d in A ])#d in days#for some reason weght 0 is not allowed and produce an spasmodic result
            date_range=[str(model.start_date+datetime.timedelta(days=d)) for d in time_range]

        #lwO heat map and point map
        lwO_values=np.moveaxis(np.array(lwO_values),0,1).tolist()
        utils.plotHeatMap(lwO_values,date_range)
        utils.plotAnimated3D(lwO_values,date_range)
        #utils.plotPointMap(lwO_values,date_range)
        #Adults map and point map
        A_values=np.moveaxis(np.array(A_values),0,1).tolist()
        utils.plotHeatMap(A_values,date_range)
        utils.plotAnimated3D(A_values,date_range)
        #utils.plotPointMap(A_values,date_range)


def exportCSV(data):
    for key in data:
        #header
        csv='Date,'
        for serie in data[key]: csv+=serie['name']+','
        csv+='\n'
        #body
        y=[ [str(datetime.date()) for datetime in data[key][0]['x']] ]
        for serie in data[key]: y+=[serie['y']]
        y=np.array(y).transpose()
        for i in range(len(y)): csv+=','.join(y[i,:]) + '\n'
        #write csv
        open('out/'+key+'.csv','w').write(csv)

try:
    from otero_precipitation_wrapper import Model as _Model
except ImportError:
    pass
import time
import tempfile
def runCpp():
    for i in range(0,100):
        configuration=Configuration('resources/otero_precipitation.cfg')
        n=configuration.getArray('breeding_site','distribution').shape[0]
        vr=np.random.rand(n)
        configuration.config_parser.set('breeding_site','distribution',', '.join( map(str,vr/vr.sum()) ))
        configuration.config_parser.set('breeding_site','height', ', '.join( map(str,vr*300 + 1) ) )
        configuration.config_parser.set('breeding_site','bare', ', '.join( map(str,vr) ) )
        configuration.config_parser.set('breeding_site','evaporation_factor', ', '.join( map(str,vr*2) ) )
        config_filename=tempfile.NamedTemporaryFile(suffix='.cfg').name
        with open(config_filename, 'w') as configfile:
            configuration.config_parser.write(configfile)

        start = time.process_time()
        model=_Model(config_filename)
        time_range,Y1=model.solveEquations()
        Y1=np.array(Y1)
        elapsed1=time.process_time() - start

        start = time.process_time()
        model=Model(Configuration(config_filename))
        time_range,Y2=model.solveEquations()
        elapsed2=time.process_time() - start

        print('||Y2-Y1||/||Y2|| * 100=%s , t2/t1=%s, %s'%(np.linalg.norm(Y2-Y1)/np.linalg.norm(Y2) *100, elapsed2/elapsed1,config_filename) )

#import netCDF4 as nc
def runInfo(nc_filename):
    grp = nc.Dataset(nc_filename)
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    lat=-31.420082
    lon=-64.188774
    p=precipitations[(abs(lons-lon)).argmin(),(abs(lats-lat)).argmin()]
    print(precipitations.shape)
    print(grp)
    print(p)

def rewritehistory():
    DATA_FOLDER='data/public/'
    HISTORY_FOLDER=DATA_FOLDER  + '.history/'
    for filename in  os.listdir(HISTORY_FOLDER):
        if(filename=='.empty' or not open(HISTORY_FOLDER+filename).readline().startswith('Date')): continue
        location,year,month,day=filename.replace('.full.weather','').replace('.csv','').split('-')
        start_date,tmp=utils.getStartEndDates(HISTORY_FOLDER+filename)
        end_date=datetime.date(int(year),int(month),int(day))
        precipitations=utils.getPrecipitationsFromCsv(DATA_FOLDER+location+'.full.csv',start_date,end_date)
        content='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h\n'
        for i,line in enumerate(open(HISTORY_FOLDER+filename).readlines()[1:]):
            fields=line.rstrip().split(',')
            if(i<len(precipitations)): fields[4]=str(precipitations[i])
            content+=','.join(fields)+',,\n'
        open('data/public/out/'+filename,'w').write(content)
        #print(start_date,end_date)

def weeks():
    date=datetime.date(2019,1,1)
    for i in range(0,54):
        print(i+1,' :',date+ datetime.timedelta(i*7))

import math
def rates():
    temperatures=np.array([22.13,20.73,18.89,12.98])+273.15

    for T_t in temperatures:
        elr,lpr,par,ovr1,ovr2=vR_D(T_t)
        print('*'*30 + str(round(T_t-273.15,2)) + '*'*30 )
        print('Development Time (in days)')
        #print('-'*65)
        print('\t Egg->Larva: %s'%(1/elr))
        print('\t Larvae->Pupa: %s'%(1/lpr))
        print('\t Pupa->Adult1: %s'%(1/par))
        print('\t Adult1->Adult2: %s'%(1/ovr1))
        print('Survival (in days)')
        #print('-'*65)
        me=0.01#mortality of the egg, for T in [278,303]
        ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
        mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
        print('\t Egg: %s'%(1/me))
        print('\t Larva: %s'%(1/ml))
        print('\t Pupa: %s'%(1/mp))
        #print('-'*65)
        print('Proportion:')
        print('\t Female: 50.0')
        print('\t Male: 50.0')
        #print('-'*65)
        print('\n')

def call_fitter(ovitrap_id):
    return os.system('python src/equation_fitter.py %s %s'%(ovitrap_id,sys.argv[2]) )

def fit():
    mp.Pool(mp.cpu_count()-4).map(call_fitter, range(1,152))

from equation_fitter import OVI_FIT
def plotFittedOvitrap(ovitrap_id,title=''):
    model=Model(utils.getFittedConfiguration( OVI_FIT%(sys.argv[2],ovitrap_id) ))
    time_range,Y=model.solveEquations()
    utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list([ovitrap_id]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
    title=title,
    xaxis_title='Date',
    yaxis_title='Eggs')

def plotFittedResults():
    lwO_values=[]
    A_values=[]
    O_PLUS_STD=260.
    A_MAX=35./15.
    best_ovi={'id':0,'fun':500}
    worst_ovi={'id':0,'fun':0}
    for ovitrap_id in range(1,152):
        if(not os.path.isfile( OVI_FIT%(sys.argv[2],ovitrap_id) )):
            print('%s, '%(ovitrap_id),end='')
            continue
        #solve the model fot the fitted parameters
        configuration=utils.getFittedConfiguration( OVI_FIT%(sys.argv[2],ovitrap_id) )
        configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
        model=Model(configuration)
        time_range,Y=model.solveEquations()
        BS_a,vBS_d,m,n,OVIPOSITION,ADULT1,ADULT2=model.parameters.BS_a,model.parameters.vBS_d,model.parameters.m,model.parameters.n,model.parameters.OVIPOSITION,model.parameters.ADULT1,model.parameters.ADULT2
        Y=model.Y
        indexOf=lambda t: (np.abs(time_range-t)).argmin()
        lwO=np.array([ (Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION]).reshape(m,n).sum(axis=0) for t in time_range])/(BS_a*vBS_d)#here index and day number are equal
        lwO=lwO[:,0]#if multiple container w assume the first one is ovitrap and the rest are wild containers
        A=(Y[:,ADULT1]+Y[:,ADULT2])/BS_a
        lat,lon=utils.getCoord('data/private/coord.csv',ovitrap_id)
        lwO_values.append([[lat,lon,lwO_d/O_PLUS_STD + 1e-4] for lwO_d in lwO ])#d in days#for some reason weght 0 is not allowed and produce an spasmodic result
        A_values.append([[lat,lon,A_d/A_MAX + 1e-4] for A_d in A ])#d in days#for some reason weght 0 is not allowed and produce an spasmodic result
        date_range=[str(model.start_date+datetime.timedelta(days=d)) for d in time_range]

        #Saver worst/best performing ovi to plot later
        fun=float(utils.extractPattern(utils.FITTED_FUN, OVI_FIT%(sys.argv[2],ovitrap_id)))
        if(fun<best_ovi['fun']):best_ovi={'id':ovitrap_id,'fun':fun}
        if(fun>worst_ovi['fun']):worst_ovi={'id':ovitrap_id,'fun':fun}

    #lwO heat map and point map
    lwO_values=np.moveaxis(np.array(lwO_values),0,1).tolist()
    utils.plotHeatMap(lwO_values,date_range)
    utils.plotPointMap(lwO_values,date_range)
    #Adults map and point map
    A_values=np.moveaxis(np.array(A_values),0,1).tolist()
    utils.plotHeatMap(A_values,date_range)
    utils.plotPointMap(A_values,date_range)
    #plot results
    plotFittedOvitrap(best_ovi['id'],best_ovi['fun'])
    plotFittedOvitrap(worst_ovi['id'],worst_ovi['fun'])
    plotFittedOvitrap(151)

if(__name__ == '__main__'):
    if(len(sys.argv)>1 and sys.argv[1]=='spatial'):
        runSpatial()
    elif(len(sys.argv)>1 and sys.argv[1]=='cpp'):
        runCpp()
    elif(len(sys.argv)>1 and sys.argv[1]=='info'):
        runInfo(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='rewrite'):
        rewritehistory()
    elif(len(sys.argv)>1 and sys.argv[1]=='weeks'):
        weeks()
    elif(len(sys.argv)>1 and sys.argv[1]=='rates'):
        rates()
    elif(len(sys.argv)>1 and sys.argv[1]=='fit'):
        fit()
    elif(len(sys.argv)>1 and sys.argv[1]=='plotFit'):
        plotFittedResults()
    elif(len(sys.argv)>1 and sys.argv[1]=='plotFitConf'):
        utils.kmeansFittedConfiguration([ OVI_FIT%(sys.argv[2],ovitrap_id) for ovitrap_id in range(1,152)],clusters=3)
    else:#the default is just a number indicating which test case to run, or none (test case 1 will will be default)
        if(len(sys.argv)<2):
            case=1
        else:
            case=int(sys.argv[1])
        runCases(case)
