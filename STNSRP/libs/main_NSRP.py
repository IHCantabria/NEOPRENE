
from MathematicalPropertiesSTNSRP import *
from libs_others import *
from libs_NSRP import *
import time

def main_NSRP(Datos_, Seasonality, temporal_resolution, year_ini, year_fin, process, Intensity_function,\
             statistics, weights, number_iterations, abejas, TiempoEntreTormentas, NumeroCeldasTormenta, \
             DuracionCelda, IntensidaCelda, DesplazamientoCeldasInicioTromenta):
    

    I_F=Intensity_function
    ##Saco los límites a la resolución temporal indicada
    if temporal_resolution=='d':
        t=24
    else:
        t=1
    lim = np.array([
        [(1/TiempoEntreTormentas[1])*t, (1/TiempoEntreTormentas[0])*t],#tiempo entre tormentas
        [NumeroCeldasTormenta[0], NumeroCeldasTormenta[1]],
        [(1/DuracionCelda[1])*t, (1/DuracionCelda[0])*t],
        [(1/IntensidaCelda[1])*t, (1/IntensidaCelda[0])*t], 
        [(1/DesplazamientoCeldasInicioTromenta[1])*t, (1/DesplazamientoCeldasInicioTromenta[0])*t]])
    
    ##Comienzo con la simulación y guardo lo parámetros ajustados
    #%%capture
    ##Creo dataframe para guardar los estadísticos reales
    statististics_dataframe=pd.DataFrame(index=statistics)
    ##Estadisticos ajustados
    statististics_fit_dataframe=pd.DataFrame(index=statistics)
    ##Creo dataframe para guardar los parámetros ajustados
    statististics_simulacion_dataframe=pd.DataFrame(index=statistics)
    ##Creo dataframe para guardar los parámetros ajustados
    if process=='normal'and Intensity_function=='E':
        param_s=['landa', 'ipsilon', 'eta', 'xi', 'betha'] 
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=5
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]]])
    elif process=='normal'and Intensity_function=='W':
        param_s=['landa', 'ipsilon', 'eta', 'xi', 'betha', 'alpha']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=6
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 3]])  
    elif process=='normal'and Intensity_function=='G':
        param_s=['landa', 'ipsilon', 'eta', 'xi', 'betha', 'alpha']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=6
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 3]])  
    elif process=='cells' and Intensity_function=='E':
        param_s=['landa', 'ipsilon', 'eta1', 'xi1', 'betha', 'alpha_p1', 'eta2', 'xi2', 'alpha_p2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=8
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [0, 1], [lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]]]) 
    elif process=='cells' and Intensity_function=='W':
        param_s=['landa', 'ipsilon', 'eta1', 'xi1', 'betha', 'alpha1', 'alpha_p1', 'eta2', 'xi2','alpha2', 'alpha_p2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=9
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 5], [0, 1], [lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]]]) 
    elif process=='cells' and Intensity_function=='G':
        param_s=['landa', 'ipsilon', 'eta1', 'xi1', 'alpha1', 'betha', 'alpha_p1', 'eta2', 'xi2','alpha2', 'alpha_p2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=9
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 5], [0, 1], [lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]]]) 
    elif process=='storms' and Intensity_function=='E':
        param_s=['landa1', 'ipsilon1', 'eta1', 'xi1', 'betha1', 'landa2', 'ipsilon2', 'eta2', 'xi2', 'betha2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=10
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],\
                          [lim[3,0],lim[3,1]],[lim[4,0],lim[4,1]]]) 
        #limites=[limites_[0], limites_[1], limites_[2], limites_[3], limites_[4], limites_[0], limites_[1], limites_[2], limites_[3], limites_[4]]
    elif process=='storms' and Intensity_function=='W':
        param_s=['landa1', 'ipsilon1', 'eta1', 'xi1', 'betha1', 'alpha1', 'landa2', 'ipsilon2', 'eta2', 'xi2', 'betha2', 'alpha2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=12
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 5], [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],\
                          [lim[3,0],lim[3,1]],[lim[4,0],lim[4,1]], [1,5]]) 
    elif process=='storms' and Intensity_function=='G':
        param_s=['landa1', 'ipsilon1', 'eta1', 'xi1', 'betha1', 'alpha1', 'landa2', 'ipsilon2', 'eta2', 'xi2', 'betha2', 'alpha2']
        Dataframe_parametros=pd.DataFrame(index=param_s); n_p=12
        limites=np.array([[lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],[lim[3,0],lim[3,1]],\
                          [lim[4,0],lim[4,1]], [1, 5], [lim[0,0],lim[0,1]],[lim[1,0],lim[1,1]],[lim[2,0],lim[2,1]],\
                          [lim[3,0],lim[3,1]],[lim[4,0],lim[4,1]], [1,5]]) 

    ##Comienzo el ajuste para cada con estacionalidad
    print('')
    print('')
    print('##################################################################')
    print('Adjustment of parameters using the PSO SWARM optimization function')
    print('')
    print('')
    for pri, prii in enumerate(Seasonality):
        ###Preproceso los datos
        Datos=Datos_.copy(); Datos[Datos['Rain']<0]=np.nan
        ##Estacionalid
        if len(Seasonality)==12:
            ##Selecciono solo las fechas (meses) en los que voy a calcular los estadísticos para ajustar mas tarde.
            Datos['Estacionalidad']=Datos['Rain']*np.nan
            pos=np.where(Datos.index.month == prii); pos=pos[0]
            Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]###Cojo los datos que no tienen nans
            Datos=Pluvio_GS
            #Datos.plot(figsize=(20, 8), xlim=("2005", "2006"),  ylim=(0, 80))
        else:
            ##Selecciono solo las fechas (meses) en los que voy a calcular los estadísticos para ajustar mas tarde.
            Datos['Estacionalidad']=Datos['Rain']*np.nan
            for i, ii in enumerate(prii):
                pos=np.where(Datos.index.month == ii); pos=pos[0]
                Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]###Cojo los datos que no tienen nans
            Datos=Pluvio_GS
            #Datos.plot(figsize=(20, 8), xlim=("2005", "2006"),  ylim=(0, 80))

        ##Calculo estadísticos que voy a ajustar y los meto en un dataframe
        statististics_values_real=calculate_statistics(Datos,statistics, temporal_resolution)
        statististics_dataframe[prii]=statististics_values_real
        ##Minimizo los errores en los estadísticos que quiero ajustar
        evaluator = evaluateInd_PSO(statististics_dataframe[prii].values, weights, I_F, process, statistics)
        
        print('')
        print('Fitting the months = ' + str(prii))
        print('')
        Steps = 0
        mySwarm = pso.Swarm(n_p, limites, abejas, evaluator, 'min') ##1000 abejas que se mueven en 5 direcciones para 
        ##tratar de minimizar el error de la ecuacion 27
        while True:
            print('Iteration number  = ' + str(Steps))
            mySwarm.step()
            Steps += 1
            (b, bp) = mySwarm.getBest()
            #print(evaluator.compute(bp))
            #print(bp)
            #print('Tiempo entre tromentas                     = ' + str(t*(1/bp[0])) + ' horas')
            #print('Numero de celdas por tormenta              = ' + str(bp[1]))
            #print('Duracion Celdas                            = ' + str(t*(1/bp[2])) + ' horas')
            #print('Intensidad Celdas                          = ' + str(t*(1/bp[3])) + ' mm/hora')
            #print('Desplazamiento de las celdas a la tormenta = ' + str(t*(1/bp[4])) + ' horas')    
            print('Total error = ' + str(evaluator.totalError(bp)))

            #if (evaluator.totalError(bp) < 5e-3): break
            if (evaluator.totalError(bp) < 0.001): break
            elif Steps > number_iterations: ###1000
                #print('ERROR REQUIREMENT NOT MET AFTER ' + str(number_iterations) + ' ITERATIONS')
                break
        sta_fit=evaluator.compute(bp);
        statististics_fit_dataframe[prii]= sta_fit[0]
        Fitted_parameters=mySwarm.getBest()[1]

        

        if process=='normal' and Intensity_function=='E':

            param_v=list(); 
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            Dataframe_parametros[prii]=param_v

        elif process=='normal' and Intensity_function=='W':

            param_v=list(); 
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            param_v.append(Fitted_parameters[5]);
            Dataframe_parametros[prii]=param_v

        elif process=='normal' and Intensity_function=='G':

            param_v=list(); 
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            param_v.append(Fitted_parameters[5]);
            Dataframe_parametros[prii]=param_v
            

        elif process=='storms'and Intensity_function=='E':

            param_v=list();
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            param_v.append(Fitted_parameters[5]); 
            param_v.append(Fitted_parameters[6]); 
            param_v.append(Fitted_parameters[7]); 
            param_v.append(Fitted_parameters[8]); 
            param_v.append(Fitted_parameters[9]); 
            Dataframe_parametros[prii]=param_v

        elif process=='storms'and Intensity_function=='W':

            param_v=list();
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            param_v.append(Fitted_parameters[5]); 
            param_v.append(Fitted_parameters[6]); 
            param_v.append(Fitted_parameters[7]); 
            param_v.append(Fitted_parameters[8]); 
            param_v.append(Fitted_parameters[9]); 
            param_v.append(Fitted_parameters[10]); 
            param_v.append(Fitted_parameters[11]); 
            Dataframe_parametros[prii]=param_v

        elif process=='storms'and Intensity_function=='G':

            param_v=list();
            param_v.append(Fitted_parameters[0]);
            param_v.append(Fitted_parameters[1]); 
            param_v.append(Fitted_parameters[2]); 
            param_v.append(Fitted_parameters[3]); 
            param_v.append(Fitted_parameters[4]); 
            param_v.append(Fitted_parameters[5]); 
            param_v.append(Fitted_parameters[6]); 
            param_v.append(Fitted_parameters[7]); 
            param_v.append(Fitted_parameters[8]); 
            param_v.append(Fitted_parameters[9]); 
            param_v.append(Fitted_parameters[10]); 
            param_v.append(Fitted_parameters[11]); 
            Dataframe_parametros[prii]=param_v
            
    Dataframe_parametros=allmonths(Dataframe_parametros)
    print('')
    print('')
    print('##################################################################')
    print('Synthetic simulation')
    print('')
    print('')
    ##Comienzo con la simulacion sintética y calculo estadísticos
    if process=='normal':
        [Dataframe_simulacion_unido_hour_, Dataframe_simulacion_unido_day_, Intensidad_cellss, Duracion_horas_cellss]=\
        NSRP_simulation1(Dataframe_parametros, year_ini, year_fin, temporal_resolution, I_F, process, Seasonality)
        print('Lluvia total acumulada parametros = ' + str(np.sum(Intensidad_cellss*Duracion_horas_cellss)))
        print('Lluvia total acumulada simulacion = ' + str(np.sum(Dataframe_simulacion_unido_day_)))

    elif process=='storms':

        Dataframe_parametros1=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'xi', 'betha'])
        for i in Dataframe_parametros.columns:
            Dataframe_parametros1[i]=Dataframe_parametros[i][0:5].values

        Dataframe_parametros2=pd.DataFrame(index=['landa', 'ipsilon', 'eta', 'xi', 'betha'])
        for i in Dataframe_parametros.columns:
            Dataframe_parametros2[i]=Dataframe_parametros[i][5:].values


        [Dataframe_simulacion_unido_hour_1, Dataframe_simulacion_unido_day_1, Intensidad_cellss1, Duracion_horas_cellss1]=\
        NSRP_simulation1(Dataframe_parametros1, year_ini, year_fin, temporal_resolution, I_F, 'normal', Seasonality)

        [Dataframe_simulacion_unido_hour_2, Dataframe_simulacion_unido_day_2, Intensidad_cellss2, Duracion_horas_cellss2]=\
        NSRP_simulation1(Dataframe_parametros2, year_ini, year_fin, temporal_resolution, I_F, 'normal', Seasonality)

        ##Combino los dos processos
        Dataframe_simulacion_unido_hour_=pd.DataFrame(index=Dataframe_simulacion_unido_hour_1.index)
        Dataframe_simulacion_unido_hour_['Rain']=Dataframe_simulacion_unido_hour_1.values+Dataframe_simulacion_unido_hour_2.values
        Dataframe_simulacion_unido_day_=pd.DataFrame(index=Dataframe_simulacion_unido_day_2.index)
        Dataframe_simulacion_unido_day_['Rain']=Dataframe_simulacion_unido_day_1.values+Dataframe_simulacion_unido_day_2.values

        print('Lluvia total acumulada parametros storm 1 = ' + str(np.sum(Intensidad_cellss1*Duracion_horas_cellss1))) 
        print('Lluvia total acumulada parametros storm 2 = ' + str(np.sum(Intensidad_cellss2*Duracion_horas_cellss2)))  

        print('Lluvia total acumulada parametros = ' + \
              str(np.sum(Intensidad_cellss1*Duracion_horas_cellss1)+np.sum(Intensidad_cellss2*Duracion_horas_cellss2)))
        print('Lluvia total acumulada simulacion = ' + str(np.sum(Dataframe_simulacion_unido_day_)))    
    
    
    statististics_fit_dataframe=allmonths(statististics_fit_dataframe)
    print('')
    print('')
    print('##################################################################')
    print('Results Validation')
    print('')
    print('')
    ##Caculo estadísticos realies, ajustados y simulados.
    for pri, prii in enumerate(Seasonality):
        ###Preproceso los datos
        if temporal_resolution=='d':
            Datos=Dataframe_simulacion_unido_day_.copy(); Datos[Datos['Rain']<0]=np.nan
        elif temporal_resolution=='h':
            Datos=Dataframe_simulacion_unido_hour_.copy(); Datos[Datos['Rain']<0]=np.nan
            
        ##Estacionalidad
        if len(Seasonality)==12:
            ##Selecciono solo las fechas (meses) en los que voy a calcular los estadísticos para ajustar mas tarde.
            Datos['Estacionalidad']=Datos['Rain']*np.nan
            pos=np.where(Datos.index.month == prii); pos=pos[0]
            Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]###Cojo los datos que no tienen nans
            Pluvio_GS[Pluvio_GS<0.001]=0
            Datos=Pluvio_GS
            #Datos.plot(figsize=(20, 8), xlim=("2005", "2006"),  ylim=(0, 80))
        else:
            ##Selecciono solo las fechas (meses) en los que voy a calcular los estadísticos para ajustar mas tarde.
            Datos['Estacionalidad']=Datos['Rain']*np.nan
            for i, ii in enumerate(prii):
                pos=np.where(Datos.index.month == ii); pos=pos[0]
                Datos['Estacionalidad'][pos]=Datos['Rain'][pos]
            Pluvio_GS = Datos['Estacionalidad'][Datos['Estacionalidad']>=0]###Cojo los datos que no tienen nans
            Pluvio_GS[Pluvio_GS<0.001]=0
            Datos=Pluvio_GS
            #Datos.plot(figsize=(20, 8), xlim=("2005", "2006"),  ylim=(0, 80))

        ##Calculo estadísticos que voy a ajustar y los meto en un dataframe
        statististics_values_sintetico=calculate_statistics(Datos,statistics,temporal_resolution)
        statististics_simulacion_dataframe[prii]=statististics_values_sintetico  

    if temporal_resolution=='d':
        Datos=Dataframe_simulacion_unido_day_.copy(); Datos[Datos['Rain']<0]=np.nan
    elif temporal_resolution=='h':
        Datos=Dataframe_simulacion_unido_hour_.copy(); Datos[Datos['Rain']<0]=np.nan
        
    print('Rainy days % (Observed) ' + str(100*(np.sum(Datos_>0.1)/len(Datos_))[0]))
    print('Rainy days % (Simulated) ' + str(100*(np.sum(Datos>0.1)/len(Datos))[0]))

    print('Maximun value (Observed) ' + str(Datos_.max()[0]))
    print('Maximun value (Simulated) ' + str(Datos.max()[0]))
    
    statististics_dataframe=allmonths(statististics_dataframe)
    statististics_fit_dataframe=allmonths(statististics_fit_dataframe)
    statististics_simulacion_dataframe=allmonths(statististics_simulacion_dataframe)        
			
    ##Veo si existe estacionalidad en los datos
    MK=int(ceil(sqrt(len(statistics)))); NK=int(ceil(sqrt(len(statistics))))
    fig, axes = plt.subplots(MK, NK, figsize=(22, 22))
    for i, ii in enumerate(statistics):
        Data_sta=pd.DataFrame(index=np.arange(1, 13))
        Real=list(); Ajustado=list(); Simulado=list();
        for m in np.arange(1, 13):
            Real.append(statististics_dataframe[m].loc[ii])
            Ajustado.append((statististics_fit_dataframe[m][ii]))
            Simulado.append(statististics_simulacion_dataframe[m].loc[ii])

        #ax=Dataframe_simulacion_unido_hour_1.plot(figsize=(20, 8), xlim=("2000", "2001"))
        #Dataframe_simulacion_unido_hour_2.plot(figsize=(20, 8), xlim=("2000", "2001"), style='--', color='r', ax=ax)
        Data_sta['Real']=Real
        Data_sta['Ajustado']=Ajustado
        Data_sta['Simulado']=Simulado
        #Ploteo
        col = i % MK
        row = i // MK
        ax=axes[row,col]
        ax.set_title(ii)
        #Data_sta.plot(ax=axes[row,col])
        legnd=['Observed', 'Fitted', 'Simulated']
        pp=Data_sta['Real'].plot(style='r-', lw=3,  ax=axes[row,col])
        Data_sta['Ajustado'].plot(style='bs', ax=axes[row,col])
        Data_sta['Simulado'].plot(style='g^', ax=axes[row,col])
        pp.legend(legnd, loc='best', fontsize=15)
        if ((ii[0:2]==('fi')) or (ii[0:2]==('au'))):
            ax.set_ylim(0, 1)
        #axes[row, col].text(.05, .1, str(format(i, '02d')), ha="center", va="center",
        #                  fontdict=dict(fontsize=30, fontweight='bold', color="black"), 
        #                  bbox=dict(facecolor="white", alpha=0.85), 
        #                  transform=axes[row, col].transAxes)


    #time.sleep(10) 
    ##Calculo periodos de retorno
    print('')
    print('')
    print('##################################################################')
    #print('Extreme analysis Observed')
    #RETURN_PERIODS(Datos_,  280, 'GEV', 'delta')
    #print(RETURN_PERIODS.return_periods)
    #print(RETURN_PERIODS.return_values)
    print('')
    print('')
    #print('Extreme analysis Simulated')
    #RETURN_PERIODS(Datos, 280, 'GEV', 'delta')
    #print(RETURN_PERIODS.return_periods)
    #print(RETURN_PERIODS.return_values)


    main_NSRP.Daily_Simulation=Dataframe_simulacion_unido_day_
    main_NSRP.Hourly_Simulation=Dataframe_simulacion_unido_hour_
    main_NSRP.Fitted_parameters=Dataframe_parametros
    main_NSRP.Statististics_Observed=statististics_dataframe
    main_NSRP.statististics_Fit=statististics_fit_dataframe
    main_NSRP.statististics_Simulated=statististics_simulacion_dataframe
