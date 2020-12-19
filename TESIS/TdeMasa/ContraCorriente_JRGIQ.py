# -*- coding: utf-8 -*-
"""
Created on Sat May 30 00:18:52 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.
"""

#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.
import math as mt
import scipy.optimize as op
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Longitud= 3.1 ##  Longitud de la torre [m].
DiametroTorre= 1  ##  Diámetro de la torre [m].
AreaF=mt.pi*(DiametroTorre/2)**2  ##  Área de flujo total de la torre [m^2].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.

g= 9.8  ## Constante gravitacional [m/s^2].
MMT= 92.138  ##  Masa molar del Toluneo [kmol/kg].
MMA= 18  ##  Masa molar del agua [kmol/kg].
MMN= 162.232  ##  Masa molar de la nicotina [kmol/kg].
Ro_Tolueno= 857.7  ## Densidad del tolueno puro [kg/m^3].
Ro_Agua= 995.7  ## Densidad del agua pura [kg/m^3].
Ro_Nicot= 1010  ## Densidad de la nicotina pura [kg/m^3].
TenSuperTol= 0.0291  ##  Tensión superficial de la fase dispersa TOLUENO [N/m]. 
ViscTol= 5.24825e-4  ##  Viscocidad del tolueno puro [Pa-s].
ViscAgua= 8e-4  ##  Viscocidad del agua pura [Pa-s].
ViscNicot= 3.442e-3  ##  Viscocidad de la nicotina pura [Pa-s].
Dif_NicAgua= 6.8753e-10  ##  Difusividad de la nicotina en agua [m^2/s].
Dif_NicTol= 1.62033e-9  ##  Difusividad de la nicotina en Tolueno [m^2/s].
DiamGota= 0.001  ##  Diametro de cada gota [m].
AreaGota= mt.pi*DiamGota**2
VolGota=((mt.pi*(DiamGota**3))/(6))  ##  Volumen de cada gota [m^3].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES.(Estado estacionario).(ini=Inicial, in=Entrada, out=salida).

FlujoMolR= 0.2  ##  Flujo molar de la fase refinado [kmol/s].
FlujoMolE= 0.0158  ##  Flujo molar de la fase extracto [kmol/s].

YA_in=  0.005  ##  Fracción molar de entrada del soluto en el extracto respecto al solvente libre de soluto [kmol Nicotina/ kmol Extracto].
YA_out=0.1  ##  Fracción molar de salida de soluto para el extracto.
XA_in= 0.0087  ##  Fracción molar de entrada del soluto en el refinado respecto al solvente libre de soluto [kmol Nicotina/ kmol Refinado].
XA_out= 0.0012  ##  Fracción molar de salida del soluto en el refinado respecto al solvente libre de soluto a la salida de la torre [kmol Nicotina/ kmol Refinado].



#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  ESPACIO DE SIMULACIÓN, PASO Y NÚMERO DE PARTICIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 3.1- TIEMPO DE SIMULACIÓN Y NÚMERO DE PARTICIONES EN EL TIEMPO.

Long_Inicial= 0  ##  Longitud inicial de la simulación [s].
Long_Final = Longitud  ##  Longitud final de la torre [m].
PasoEsp = 0.05 ## Longitud para cada tramo de la torre [m].
NumPasosEsp=int((Long_Final-Long_Inicial)/(PasoEsp))  ##  Determina el número de particiones necesarias para resolver el método.
Long = np.linspace(Long_Inicial, Long_Final, NumPasosEsp)  ##  Crea vector de longitud para poder graficarlo [s].
DiscCDE=0.012/NumPasosEsp  ##  Discretización de la fracción molar de soluto en el refinado a la entrada para graficar la CDE.
VolPart= AreaF*PasoEsp  ##  Valor de cada partición de volumen [m^3].
Iter= 100  ##  Número de iteraciones para lograr convergencia.
Tolerancia=1e-8  ##  Tolerancia aceptada para las iteraciones.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final = 1000#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
PasoTemp = 10 ##float(input("ingrese el tamaño del PasoTemp para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
NumPasosTemp=int((t_Final-t_Inicial)/(PasoTemp))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, NumPasosTemp)  ##  Crea vector de tiempo para poder graficarlo [s].
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

XA_LdeO= np.zeros(len(Long)+1)  ##  Crea un vector de datos para la concentración del soluto A en el refinado para la LdeO.
YA_LdeO=np.zeros(len(Long)+1)  ##  Crea un vector de datos para la concentración del soluto A en el extracto para la LdeO.

XA=np.zeros(len(Long)+1) ##  Crea un vector de datos para la concentración del soluto A en el refinado para las concentraciones en equilibrio.
YA=np.zeros(len(Long)+1)  ##  Crea un vector de datos para la concentración del soluto A en el extracto para las concentraciones en equilibrio.

XA_Sup=np.zeros(len(Long))  ##  Crea un vector de datos para la concentración del soluto A para en el refinado, esta es la supuesta para poder hallar la verdadera concentración en el equilibrio.
XA_Part=np.zeros(len(Long)+1)  ## ##  Crea un vector de datos para la partición de la concentración del soluto A en el refinado.


AM=np.zeros(len(Long))  ## Crea un vector de datos para calcular el área de transferencia de masa.
NA_EE=np.zeros(len(Long))  ## Crea un vector de datos para calcular el flux de transferencia de masa.

XA_Equ=np.zeros(len(Long))  ##  Crea un vector de datos para el calculo de la concentración del soluto A (refinado) en equilibrio o interfaz.
YA_Equ=np.zeros(len(Long))  ##  Crea un vector de datos para el calculo de la concentración del soluto A (extracto) en equilibrio o interfaz.

XA_CDE=np.zeros(len(Long))  ##  Crea un vector de datos para la concentración del soluto A en la fase refinado para graficar la CDE.
YA_CDE=np.zeros(len(Long))  ##  Crea un vector de datos para la concentración del soluto A en la fase extracto para graficar la CDE.


#________________________________________________________________________________________________________________________________________________________________________
## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA.##


XA[0]=XA_in  ##  Otorga un valor inicial en la posición cero del vector para la concentracion de A en la fase refinado.
YA[0]=YA_out  ##  Otorga un valor inicial en la posición cero del vector para la concentracion de A en la fase extracto.
XA[NumPasosEsp]=XA_out  ##  Otorga el valor de salida o final en la última partición de la torre para la fase refinado.

XA_LdeO[0]= XA_in  ##  Otorga el valor inicial a la concentración de refinado para iniciar con la construcción de la linea de operación.
YA_LdeO[0]= YA[0] ##  Otorga el valor inicial a la concentración de extracto para iniciar con la construcción de la linea de operación.

mLdeO=((YA_in-YA[0])/(XA[NumPasosEsp]-XA_in))  ##  Cálcula la pendiente de la linea de operación en funcion de la variación de concentraciónes del extracto y refinado.

XA_IntTotal= XA_in-XA[NumPasosEsp]  ##  Intervalo de concentración del refinado desde la entrada hasta la salida con concentración deseada [kmol Nicotina/ kmol Refinado].

for j in range(0,NumPasosEsp,1): 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
#  MODULO 5-  SOLUCIÓN DE ECUACIONES.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
#  MODULO 5.1-  CREACIÓN DE DATOS INICIALES PARA LAS MATRICES DE CONCENTRACIÓN.
    

    
    XA_LdeO[j+1]=XA_LdeO[j]-(XA_IntTotal/NumPasosEsp)  ##  Otorga el siguiente valor a la concentración en el refinado para construir la linea de operación.
    YA_LdeO[j+1]=mLdeO*(XA_LdeO[j+1]-XA[NumPasosEsp])+YA_in  ##  Otorga el siguiente valor a la concentración en el extracto para construir la linea de operación.


    XA_Sup[j]=XA_in-j*(XA_IntTotal/NumPasosEsp)  ##  Se dan valores semilla para la concentración de A en la fase refinado. 

    if j==0:  ##  Condiciona la concentración de cada partición para que la siguiente comience con la salida de la última.
        XA_Part[j]=XA_in        
    else:
        
        XA_Part[j]=XA[j-1]
        
    if j==0:
        YA[j]=YA[0]

    
    ErrorEst= 2*Tolerancia  ## Se define un error estimado.

    Contador=0  ##  Contador para ingresar al ciclo while.

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
## CICLO WHILE ##  Este ciclo se crea para encontrar las concentraciones en equilibrio contando con valores semillas para cada partición.
    
    while ErrorEst>Tolerancia or  Contador<Iter:
        
                
  
        xA_Sup=((XA_Sup[j])/(1+XA_Sup[j]))  ##  Fracción molar respecto de la fase como un todo.
        yA=((YA[j])/(1+YA[j]))
                
        XwA_Sup=XA_Sup[j]*(MMN/MMA)  ##  Conversión a fracción molar respecto del solvente puro.
        YwA_Sup=YA[j]*(MMN/MMT)
                
        xwA_Sup=((XwA_Sup)/(1+XwA_Sup))  ##  Conversión a fracción másica respecto de la fase como un todo.
        ywA_Sup=((YwA_Sup)/(1+YwA_Sup))
                       
    
        Rho_R=1/((xwA_Sup/Ro_Nicot)+((1-xwA_Sup)/(Ro_Agua)))  ##  Cálculo y actualizacion de las densidades.
        Rho_E=1/((ywA_Sup/Ro_Nicot)+((1-ywA_Sup)/(Ro_Tolueno)))
                
        MMR_Mezcla=MMA*(1-xA_Sup)+MMN*xA_Sup ##  Cálculo y actualizacion de las masas molares.
        MME_Mezcla=MMT*(1-yA)+MMN*yA        
                
        ViscR= np.exp(xwA_Sup*np.log(ViscNicot)+(1-xwA_Sup)*np.log(ViscAgua))  ##  Cálcuo y actualización del las viscosidades.
        ViscE= np.exp(ywA_Sup*np.log(ViscNicot)+(1-ywA_Sup)*np.log(ViscTol))
                
        FMolTotR=FlujoMolR-FlujoMolR*XA_Sup[j]  ##  Cálcula y actualiza el flujo molar de la fase refinado.
        FMasTotR=FMolTotR*MMR_Mezcla  ##  Cálcula y actualiza el flujo masico de la fase refinado.
        FVolTotR=FMasTotR/Rho_R  ##  Cálcula y actualiza el flujo volumétrico de la fase refinado.
                
        FMolTotE=FlujoMolE+FlujoMolE*YA[j]  ##  Cálcula y actualiza el flujo molar de la fase extracto. 
        FMasTotE=FMolTotE*MME_Mezcla  ##  Cálcula y actualiza el flujo masico de la fase extracto.
        FVolTotE=FMasTotE/Rho_E  ##  Cálcula y actualiza el flujo volumétrico de la fase extracto.
        
        FracVolE=((FVolTotE)/(FVolTotR+FVolTotE))  ##  Cálcula y actualiza la fracción volumétrica en la fase extracto.
        VelR=((FVolTotR)/(AreaF*(1-FracVolE)))  ##  Cálcula y actualiza la velocidad de la fase refinado.
        VelE=((FVolTotE)/(AreaF*FracVolE))  ##  Cálcula y actualiza la velocidad de la fase extracto.
        VelFases=VelR+VelE  ##  Cálcula y actualiza la velocidad de deslizamiento de las fases.
                
        ReGota= ((DiamGota*VelFases*Rho_R)/(ViscR))  ##  Cálcula y actualiza el número de Reynolds de la gota.
        ScR=((ViscR)/(Rho_R*Dif_NicAgua))  ##  Cálcula y actualiza el número de Schmidt en el refinado.
        K_Skll=(ReGota**(1/8))*((ViscR/ViscE)**(1/4))*(((ViscR*VelFases)/(TenSuperTol))**(1/6))  ##  Cálcula y actualiza la correlación de A.H.P. Skelland.
        F_Skll=0.281+(1.615*K_Skll)+(3.73*K_Skll**2)-(1.874*K_Skll**3)  ##  Cálcula y actualiza la corrección de forma esférica de la correlacion de A.H.P. Skelland.
        kX= (Dif_NicTol/DiamGota)*(2+0.463*(ReGota**0.484)*(ScR**0.339)*(((DiamGota*(g**(1/3)))/(Dif_NicAgua**(2/3)))**(0.072)))*F_Skll  ##  Cálcula y actualiza el coeficiente local de TdeMasa en la fase refinado.
        kY=(0.00375*VelFases/(1+(ViscE/ViscR)))  ##  Cálcula y actualiza el coeficiente local de TdeMasa en la fase extracto.
        
        mFI=-kY/kX  ##  Cálcula y actualiza la pendiente de la fuerza impulsora (FI).
            
        
            
        f = lambda x: -YA[j]-mFI*(XA_Sup[j]-x)+16.07781726*x+((0.08593565238*x)/(x+0.000336023425))  ## Se define el polinomio que representa la CDE.
        
        XA_Equ[0]=XA_Sup[0]  ##  Se inicia el método de Newton para encontrar las raices del polinomio y asi hallar el valor de las concentraciones en el equilibrio.
            
        XA_Equ[j]= op.newton(f,XA_Equ[j],tol=1e-8,maxiter=100)  ##  Ejecuta el método de Newton.
        
        YA_Equ[j]=16.07781726*XA_Equ[j]+((0.08593565238*XA_Equ[j])/(XA_Equ[j]+0.000336023425))    ##  Cálcula la concentracion de la fase extracto en el equilibrio.
        
        NA_EE[j]=kX*(XA_Sup[j]-XA_Equ[j])  ##  Cálcula el flux de transferencia de masa en el equilibrio.
        
        Num_Got=FracVolE*(VolPart/VolGota)  ##  Determina el número de gotas de la fase dispersa.
        AM[j]=Num_Got*AreaGota  ##  Cálcula el área de transferencia de masa.
        XA[j]=(XA_Part[j]*FlujoMolR-AM[j]*NA_EE[j])/FlujoMolR  ##  Cálcula la concentración de A en el refinado en la línea de operación.
        YA[j+1]=(YA[j]*FlujoMolE-AM[j]*NA_EE[j])/FlujoMolE   ##  Cálcula la concentración de A en el extracto en la línea de operación.
        
            
        ErrorEst=abs(XA[j]-XA_Sup[j])  ##  Analiza el error entre la concentración supuesta y la calculada por el método.
        Contador=Contador+1  ##  Aumenta de uno en uno el contador del ciclo while.
    
        XA_Sup[j]=XA[j]  ##  Condición para que en dado caso la concentración supuesta no sea la misma cálculada, ingrese de nuevo al ciclo con la concentración calculada.
                
#--------------------------------------------------------------------------------------------------------------------#
### ECUACIONES PARA LA CDE######################

    XA_CDE[j]=j*DiscCDE  ##  Gráfica la concentración del soluto A en el refinado en la CDE.
    YA_CDE[j]=16.07781726*XA_CDE[j] + (0.08593565238*XA_CDE[j])/(XA_CDE[j]+0.000336023425) ##  Gráfica la concentración del soluto A en el extracto en la CDE.
    
    
    
pt.figure("TdeMasa contracorriente", [10,5])  ##  Crea figura y le otorga dimensiones.
pt.plot(XA,YA,'k',label='LdeO',linewidth=2)  ##  Grafica el tiempo vs Concentración deseada.
pt.xlabel("Refinado")  ##  Agrega título al eje x.
pt.ylabel('Extracto')  ##  Agrega título al eje y.
pt.title('LdeO y CDE')  ##  Agrega título en la parte superior.
pt.grid(True)  ##  Agrega la cuadrilla a la imagen
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 7})
leg.get_frame().set_alpha(0.5)  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).
pt.show()
pt.plot(XA_CDE,YA_CDE,'b',label='CDE')
pt.legend(loc='best') 
## HASTA ACÁ FUNCIONA EL PRIMER SIMULADOR, EL QUE DETERMINA LAS PRIMERAS CONCENTRACIONES EN ESTADO ESTACIONARIO.
    
#--------------------------------------------------------------------------------------------------------------------#
    
# ####### MODO DINÁMICO
    
    
    
    
    
    
#     XA_in_Dis=XA_EE[0]
#     YA_out_Dis=YA_EE[0]
    
#     XA_out_Dis=XA_EE[NumPasosEsp-1]
    
#     YA_in_Dis=YA_EE[NumPasosEsp-1]
    
#     XA_Vector[0:]=XA_EE[1:NumPasosEsp+1]
#     YA_Vector[0:]=YA_EE[0:NumPasosEsp]
    
#     XA_Ent=XA_EE[0]
#     YA_Sal=YA_EE[0]
    
#     Tolerancia=1e-8
    
#     Int_X_A=XA_in_Dis-XA_out_Dis
    
#     m_LdeO=((YA_in_Dis-YA_out_Dis)/(XA_out_Dis-XA_in_Dis))

#     X_LdeO[0]=XA_in_Dis
#     Y_LdeO[0]=YA_out_Dis
    
#     X_LdeO[j+1]=X_LdeO[j]-(Int_X_A/NumPartEsp)
#     Y_LdeO[j+1]=m_LdeO*(X_LdeO[j+1]-XA_out_Dis)+YA_in_Dis
    

# X_A[0][:]=XA_Vector
# Y_A[0][:]=YA_Vector
# A_M[0][:]=AM
# N_A[0][:]=NA_EE

# XA_In[0]=XA_in_Dis
# YA_In[0]=YA_in_Dis




# FlujoM_R=FlujoMolR+FlujoMolR*XA_in_Dis

# MolesR=(VolPart*Ro_Agua)/MMA
# MolesE=(VolPart*Ro_Tolueno)/MMT


# for i in range(0,NumPasosTemp,1):
    
#     YA_In[i+1]=YA_in_Dis
#     XA_In[i+1]=XA_in_Dis
    
    

    
#     for k in range(0,NumPasosEsp,1):
                
#         x_A[i][k]=((X_A[i][k])/(1+X_A[i][k]))
#         y_A[i][k]=((Y_A[i][k])/(1+Y_A[i][k]))
        

#         Xw_A=X_A[i][k]*(MMN/MMA)  ##  Conversión a fracción molar respecto del solvente puro.
#         Yw_A=Y_A[j][k]*(MMN/MMT)
                
#         xw_A=((Xw_A)/(1+Xw_A))  ##  Conversión a fracción másica respecto de la fase como un todo.
#         yw_A=((Yw_A)/(1+Yw_A))
                           
#         MMR_Mix=MMA*(1-x_A[i][k])+MMN*x_A[i][k]
#         MME_Mix=MMT*(1-y_A[i][k])+MMN*y_A[i][k]       
                
#         Visc_R= np.exp(x_A[i][k]*np.log(ViscNicot)+(1-x_A[i][k])*np.log(ViscAgua))
#         Visc_E= np.exp(y_A[i][k]*np.log(ViscNicot)+(1-y_A[i][k])*np.log(ViscTol))        
            
#         Ro_R=1/((xw_A/Ro_Nicot)+((1-xw_A)/(Ro_Agua)))
#         Ro_E=1/((yw_A/Ro_Nicot)+((1-yw_A)/(Ro_Tolueno)))
        
#         if k==0:
            
            
#             FlujoR[i][k]=FlujoM_R-A_M[i][k]*N_A[i][k]
            
#         else:
            
#             FlujoR[i][k]=FlujoR[i][k-1]-A_M[i][k]*N_A[i][k]
            
            
#         FlujoMas_R[i][k]=FlujoR[i][k]*MMR_Mix
#         FlujoV_R=FlujoM_R/Ro_R
        
#         FlujoE[i][k]=FlujoMolE+FlujoMolE*Y_A[i][k]
        
#         FlujoMas_E[i][k]=FlujoE[i]*MME_Mix
#         FlujoV_E=FlujoMas_E[i]/Ro_E           
        
        
#         ########### Python me arroja un error porque en la línea 350 se cálcula FlujoE[i][k]
#                                                ## con las dos dimensiones pero en la línea 353 se utiliza unicamente con la dimensión temporal.
#                                                ### Python no me deja hacer esto.
        
        
        
#         FracVol_E=((FlujoV_E)/(FlujoV_R+FlujoV_E))
#         Vel_R=((FlujoV_R)/(AreaF*(1-FracVol_E)))
#         Vel_E=((FlujoV_E)/(AreaF*FracVol_E))
#         Vel_Fases=Vel_R+Vel_E
                
#         Re_Gota= ((DiamGota*Vel_Fases*Ro_R)/(Visc_R))
#         S_cR=((Visc_R)/(Ro_R*Dif_NicAgua))
#         KSkll=(Re_Gota**(1/8))*((Visc_R/Visc_E)**(1/4))*(((Visc_R*Vel_Fases)/(TenSuperTol))**(1/6))  
#         FSkll=0.281+(1.615*KSkll)+(3.73*KSkll**2)-(1.874*KSkll**3)
#         k_X= (Dif_NicTol/DiamGota)*(2+0.463*(Re_Gota**0.484)*(S_cR**0.339)*(((DiamGota*(g**(1/3)))/(Dif_NicAgua**(2/3)))**(0.072)))*FSkll
#         k_Y=(0.00375*Vel_Fases/(1+(Visc_E/Visc_R)))
        
#         m_FI=-k_Y/k_X   
        
        
        
#         f = lambda x: -Y_A[i][k]-m_FI*(X_A[i][k]-x)+16.07781726*x+((0.08593565238*x)/(x+0.000336023425))

#         if k==0:
            
#             x0=0.001
            
#         else:
            
#             x0=XA_Equil[i][k-1]
            


            
#         XA_Equil[i][k]= op.newton(f,x0,tol=1e-8,maxiter=100)
        
#         YA_Equil[i][k]=16.07781726*XA_Equil[i][k]+((0.08593565238*XA_Equil[i][k])/(XA_Equil[i][k]+0.000336023425))  
        
#         NA_EE2[i][k]=k_X*(X_A[i][k]-XA_Equil[i][k])
        
#         Num_Got=FracVol_E*(VolPart/VolGota)
#         AM2[i][k]=Num_Got*AreaGota
        
        
#         if k==0:
            
#             DeltaXA=(XA_In[i]*FlujoMolR-A_M[i+1][k]*NA_EE2[i+1][k]-X_A[i][j]*FlujoMolR)/MolesR
            
#         else:
            
#             DeltaXA=(X_A[i][k-1]*FlujoMolR-A_M[i+1][k]*NA_EE2[i+1][k]-X_A[i][j]*FlujoMolR)/MolesR
            
#         X_A[i+1][k]=X_A[i][k]+DeltaXA*PasoTemp
        
        # print(N_A)

        # for m in range(0,NumPasosEsp,1):
            
        #     k=NumPasosEsp+1-m
            
        #     if k==NumPasosEsp:                
                
        #         DeltaYA=(YA_in_Dis*FlujoMolE+A_M[i][k]*NA_EE2[i][k]-Y_A[i][j]*FlujoMolE)/MolesE
                
        #     else:
        #         DeltaYA=(Y_A[i][k+1]*FlujoMolE+A_M[i][k]*NA_EE2[i][k]-Y_A[i][j]*FlujoMolE)/MolesE
        
        
        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
