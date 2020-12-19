# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:22:54 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.


#________________________________________________________________________________________________________________________________________________________________________


*** C O N T R O L A D O R  D E  N I V E L  Y  P R E S I Ó N  P A R A  U N A  C A L D E R A   P I R O T U B U L A R (P ,P I ,P D  Y  P I D).***

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO***


- Se supondra una caldera de geometría cilíndrica ubicada horizontalmente.

- La variable de proceso 1, es la presión de vapor dentro de la calder, esta se regula con energía Q_i.

- La variable de proceso 2, es el nivel de la calder, este se regula con el flujo de líquido que ingresa a la caldera FlujoM_Liq_in_i .

- Se tiene en cuenta que la variación en el nivel depende del radio de la caldera porque al ser un cilindro horizontal su longitud es constante.

- Los dos lazos de control comparten una variable, te es el flujo de líquido que se evapora o flujo de vapor inicial FlujoM_LiqEvap_i.

- Al perturbar el sistema de presión sacando mas flujo de vapor a la salida, se genera una acción de control que suministra calor al sistema.

- Al suministrar calor al sistema, este evapora líquido, por ende baja el nivel perturbando asi el otro sistema, el de nivel.

- Al perturbar el sistema de nivel, se genera una acción de control que suministra flujo de líquido a la caldera.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***


Fluido de proceso = Agua líquida a T=25 [°C].

-Variable de proceso (y) = Nivel de agua en la caldera [m]."Nivel_i".

-Variable de proceso (y) = Presión de vapor dentro de la caldera.P_Vap_i.

-Variable de control "Nivel" (u) = Flujo de líquido que ingresa a la caldera.FlujoM_Liq_in_i. 

-Variable de control "Presión" (u) = Calor suministrado a la caldera desde una fuente externa.Q_i

-Variable de perturbación "Nivel" (d) = Flujo de líquido evaporado o flujo de vapor inicial.FlujoM_LiqEvap_i.

-Variable de perturbación "Presión" (d) = Flujo de vapor que sale de la caldera.FlujoM_Vap_out_i.

**PARA MAS DETALLES VER MANUAL DE USUARIO**
#________________________________________________________________________________________________________________________________________________________________________
"""
#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.
import math as mt

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

H=3  ##  Longitud de la caldera [m].
Diam=1  ##  Diámetro de la caldera [m].
Rad=Diam/2  ##  Radio de la caldera [m]
Vol_Tot=H*mt.pi*(Diam/2)**2  ##  Volumen total de la caldera [m^3].

SP_N=Rad  ##  Set point para el nivel [m].
SP_P=10  ##  Set point para la presión de vapor [Bar].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARÁMETROS DE LAS Y CONSTANTES.

Dens=1000  ##  Densidad del agua líquida [kg/m^3].
R=8.314  ##  Constante de los gases [KJ/Kmol K].
Teb=100+273.15  ##  Temperatura de ebullición del agua [K].
M=18  ##  Masa molar del agua [kg/kmol].

Al=8.712  ##  Constante A para el calculo del Cp (Líquido). (Tabla C3. Smith Van Ness).
Bl=1.25*1e-03  ##  Constante B para el Cp (Líquido). (Tabla C3. Smith Van Ness).
Cl=-0.18*1e-06  ##  Constante C para el Cp (Líquido). (Tabla C3. Smith Van Ness).

Av=3.470  ##  Constante A para el calculo del Cp (Vapor). (Tabla C1. Smith Van Ness).
Bv=1.450*1e-03  ##  Constante B para el calculo del Cp (Vapor). (Tabla C1. Smith Van Ness).
Cv=0*1e-06  ##  Constante C para el calculo del Cp (Vapor). (Tabla C1. Smith Van Ness).
Dv=0.121*1e5  ##  Constante D para el calculo del Cp (Vapor). (Tabla C1. Smith Van Ness).

A = 10.01  ##  Constante A para la ecuación de Antoine. (Propierties of gases and liquids).
B = 1611  ##  Constante B para la ecuación de Antoine. (Propierties of gases and liquids).
C = -51.4  ##  Constante C para ecuación de Antoine. (Propierties of gases and liquids).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES.(Estado estacionario).(ini=Inicial, in=Entrada, out=salida).

T_in = 20+273.15  ##  Temperatura de entrada del líquido a la caldera [K].
Nivel_ini=Rad  ##  Nivel del tanque, este inicia suponiendo que la caldera es horizontal y esta llena en un 50% [m].
P_Vapor_ini=SP_P ##  Presión de vapor inicial de la caldera [bar].
FlujoM_Vap_out=1  ##  Flujo másico de vapor inicial que sale de la caldera [kg/s].

FlujoM_Liq_in=0.97849297  ##  Flujo másico de líquido inicial que ingresa a la caldera (u) para el nivel [Kg/s].
FlujoM_L_Max=2*FlujoM_Liq_in  ##  Flujo másico máximo que puede tener el líquido al ingresar a la caldera [kg/s].

Q_ini=2785.83  ##  Flujo de energía inicial (u) para la presión de vapor [kw].
Q_Max=2*Q_ini  ##  Flujo de energía máxima que puede ingresar al proceso [kw].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN, PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final =5000#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
Paso = 1 ##float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
Iter = N_Part  ##  Número de iteraciones, se puede decir que es igual al número de particiones.
t_Pert = int(0.3 * N_Part)  ##  Tiempo y valor en que se aplica la perturbación [%].
n=1#int(input("Ingrese el tiempo de muestreo del controlador"))

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

Nivel_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el nivel en la caldera, variable de proceso (y).
P_Vap_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la presion de vapor dentro de la caldera, variable de proceso (y).
SP_N_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point del nivel.
SP_P_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point de la presión de vapor.
Q_i=np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para el calor o energía suministrad, variable de control para la presión de vapor.
FlujoM_Liq_in_i=np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para el flujo de líquido (variable de control). 
FlujoM_LiqEvap_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo másico de vapor o líquido que se evapora.
FlujoM_Vap_out_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo másico de vapor que sale de la caldera, perturbación para el control de presión(d).
Vol_Liq_i=np.zeros(len(tiempo)) ## Crea un vector para el volumen que ocupa el líquido en la caldera.
Vol_Vap_i=np.zeros(len(tiempo)) ##  Crea un vector de datos para el volumen que ocupa el vapor en la caldera.
Error_N_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para guardar el error para el nivel.
Error_P_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para guardar el error para la presión.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 5-  PUNTO DE AJUSTE (SET POINT) Y PERTURBACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.1-  SET POINT.

SP_N_i[t_Inicial:]=SP_N ##  Llena el vector del Set Point con el valor asignado como referencia para el nivel.
SP_P_i[t_Inicial:]=SP_P ##  Llena el vector del Set Point con el valor asignado como referencia para la presión de vapor.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.2-  PERTURBACIONES.

FlujoM_Vap_out_i[t_Inicial:] = FlujoM_Vap_out  ##  Inicia el valor de apertura para la válvula del vapor de salida en estado estacionario [kg/s].
FlujoM_Vap_out_i[t_Pert:] = FlujoM_Vap_out*1.5 ## Cambia el valor de apertura de la válvula de vapor a la salida. Tiempo en que inicia la perturbación t_Pert).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  PARÁMETROS DEL CONTROLADOR.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 6.1-  CONTROL DE NIVEL.

Kp_N=500  ##  Parámetro de ganancia proporcional [Acción de control que se produce [m^3/s  / m].
ti_N=10 ##  Parámetro para el tiempo integral [s].
td_N=0.01  ##  Parámetro para el tiempo derivativo [s].
e_1kN=0  ## Error en el instante de tiempo anterior (i-1).
e_2k=0  ## Error en el instante de tiempo anterior (i-2).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 6.1-  CONTROL DE PRESIÓN.

Kp_P= 0.05 ##  Parámetro de ganancia proporcional [kj/Bar].
ti_P=0.01 ##  Parámetro para el tiempo integral [s].
td_P=1  ##  Parámetro para el tiempo derivativo [s].
e_1kP=0  ##  Error en el instante de tiempo anterior (i-1).
e_2kP=0  ##  Error en el instante de tiempo anterior (i-2).
ITC=1  ##  Intervalo de tiempo de la acción de control.

#________________________________________________________________________________________________________________________________________________________________________

#-------------------------------------------------------------------------## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA DEL MODELO DEL PROCESO CON SU CONTROLADOR.##

Nivel=Nivel_ini  ##  Asigna el valor inicial al nivel del tanque para iniciar el método.
P_Vapor=P_Vapor_ini  ##  Asigna el valor inicial al la presión de vapor para iniciar el método.

for i in range(0,Iter,1):
#-----------------------------------------------------------------------------------------------
#  MODULO 7-  ACTUALIZACIÓN DE LAS VARIABLES DE PROCESO, CONTROL Y ERROR.

    Nivel_i[i] = Nivel  ##  Guarda los datos para cada solución de la Ecuación diferencial de Temperatura.
    Error_N_i[i]= SP_N_i[i]-Nivel  ##  Guarda el error en cada iteración para el control de nivel.
    FlujoM_Liq_in_i[t_Inicial]=FlujoM_Liq_in ## Inicia el valor de flujo de líquido en estado estacionario (NIVEL).[kg/s].
    
    P_Vap_i[i] = P_Vapor  ##  Guarda los datos para cada solución de la Ecuación diferencial de Presión de vapor.
    Error_P_i[i]=SP_P_i[i]-P_Vapor  ##  Guarda el error en cada iteración para el control de presión de vapor.
    Q_i[t_Inicial]=Q_ini ## Inicia el valor de flujo de energía en estado estacionario [kw].(PRESIÓN).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 8-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.1-  SOLUCIÓN DE ECUACIONES CONSTITUTIVAS.

    Tsat= (((-B)/(np.log10(P_Vap_i[i])-A))-C)+273.15  ##  Soluciona y actualiza la ecuación de Antoine(1888) para la temperatura de saturación.The propierties of gases and liquids [K].
    H_Vap=-0.0006*P_Vapor**3+0.1689*P_Vapor**2-20.951*P_Vapor+2284.4  ##  Polinomio que cálcula el calor latente en función de la Presión de saturación [kj/kg].(Con datos de la tabla A12. Termodinámica de Wark).
    
    Cp_L=(Al+(Bl*Tsat)+(Cl*Tsat**2))*((R)/(M))  ##  Soluciona y actualiza la ecuación para el Cp del líquido [Kj/kg K].(Tabla C3. Smith Van Ness).
    Cp_V=(Av+(Bv*Tsat)+(Cv*Tsat**2)+(Dv*Tsat**-2))*((R)/(M))  ##  Soluciona y actualiza la ecuación para el Cp del vapor [kj/kg K].(Tabla C1. Smith Van Ness).
    
    Vol_Liq_i [i]=(((mt.pi*Rad**2)/(2))+(Nivel_i[i]-Rad)*(np.sqrt(2*Rad*Nivel_i[i]-Nivel_i[i]**2))+(mt.asin(((-Rad+Nivel_i[i])/(Rad))))*(Rad**2))*H  ##  Soluciona y actualiza el volumen que ocupa el líquido dentro de la caldera [m^3].
    Vol_Vap_i[i]=Vol_Tot-Vol_Liq_i [i]  ##  Soluciona y actualiza el volumen que ocupa el vapor dentro de la caldera [m^3].

    FlujoM_LiqEvap_i[i]=((Q_i[i])/((Cp_L*(Teb-T_in))+H_Vap+(Cp_V*(Tsat-T_in))))  ##  Soluciona y actualiza la ecuacion de flujo másico para el vapor que ingresa [kg/s].   

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2-  SOLUCIÓN DE ECUACIONES DIFERENCIALES APLICANDO EL MÉTODO EULER.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2.1-  ECUACIÓN DIFERENCIAL PARA EL NIVEL.
    
    dLdt=((FlujoM_Liq_in_i[i]-FlujoM_LiqEvap_i[i])/(Dens*H))*((1)/(np.sqrt(2*Nivel*Rad-Nivel**2)+((Nivel-Rad)*(2*Rad-2*Nivel)/(2*np.sqrt(2*Rad*Nivel-Nivel**2)))+((Rad)/(np.sqrt(1-((-Rad+Nivel)/(Rad))**2)))))  ##  Soluciona la derivada de la ED para el nivel.
    Nivel=Nivel_ini+dLdt*Paso  ##  Soluciona la ecuación diferencial, agregandole al nivel inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    Nivel_ini=Nivel  ##  Actualiza el nivel para iniciar con el siguiente paso [m].
           
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2.2-  ECUACIÓN DIFERENCIAL PARA LA PRESIÓN DE VAPOR.
    
    dPdt=(((R*Tsat)/M)*(FlujoM_LiqEvap_i[i]-FlujoM_Vap_out_i[i])+P_Vapor*(H*dLdt*(((np.sqrt(2*Rad*Nivel-Nivel**2))+(Nivel-Rad)*(2*Rad-2*Nivel))/(2*np.sqrt(2*Rad*Nivel-Nivel**2)))+((Rad)/(np.sqrt(1-(((-Rad+Nivel)/(Rad))**2))))))/(Vol_Vap_i[i]*100)  ##  Soluciona la derivada de la ED para la presión de vapor. 
    P_Vapor=P_Vapor_ini+dPdt*Paso  ##  Soluciona la ecuación diferencial, agregandole a la presión de vapor inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    P_vapor_ini=P_Vapor  ##  Actualiza la presión de vapor para iniciar con el siguiente paso [Bar].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 9-  SELECCIÓN Y ACCIONES DEL CONTROLADOR.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MODULO 9.1-  CONTROLADOR (P) PARA EL NIVEL.

    if (td_N==0 and ti_N>=9999):  ##  Control (P), depende de los parametros td y ti.
        Delta_u = Kp_N*ITC*Error_N_i[i]  ##  Ecuación para el control (P).
                
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.2-  CONTROLADOR (PI),(PD) Y (PID) PARA EL NIVEL.
    else:
        Delta_u= Kp_N*ITC*(Error_N_i[i]-e_1kN) + Kp_N*(ITC/ti_N)*e_1kN + Kp_N*(td_N/ITC)*(Error_N_i[i]-2*e_1kN+e_2k)  ##  Ecuación para el control (PID)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.3-  CONTROLADOR (P) PARA LA PRESIÓN DE VAPOR.
        
    if (td_P==0 and ti_P>=9999):  ##  Control (P), depende de los parametros td y ti.
        Delta_uP = Kp_P*ITC*Error_P_i[i]  ##  Ecuación para el control (P).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.4-  CONTROLADOR (PI),(PD) Y (PID) PARA LA PRESIÓN DE VAPOR.
    else:                
        Delta_uP= Kp_P*ITC*(Error_P_i[i]-e_1kP) + Kp_P*(ITC/ti_P)*e_1kP + Kp_P*(td_P/ITC)*(Error_P_i[i]-2*e_1kP+e_2kP)  ##  Ecuación para el control (PID)        

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.5-  ACTUALIZACION DEL ERROR.
         
    e_1kN=Error_N_i[i]  ##  Actualiza el error anterior para el nivel (i-1).
    e_2k=e_1kN  ##  Actualiza el error anterior para el nivel (i-2).
    
    e_1kP=Error_P_i[i]  ##  Actualiza el error anterior para la presión de vapor (i-1).
    e_2kP=e_1kP  ##  Actualiza el error anterior para la presión de vapor (i-2).    

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.6-  ACCIÓNES DE CONTROL.    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.6.1-  ACCIÓN DE CONTROL PARA EL NIVEL.  
    
    if (FlujoM_Liq_in_i[i]+Delta_u)>FlujoM_L_Max: 
        FlujoM_Liq_in_i[i+1]=FlujoM_L_Max  ##  Corrige la acción de control para que no pase de su límite máximo.
        
    elif (FlujoM_Liq_in_i[i]+Delta_u)<0:
        FlujoM_Liq_in_i[i+1]=0  ##  Corrige la acción de control para que no pase de su límite mínimo.
    
    else:
          FlujoM_Liq_in_i[i+1]=FlujoM_Liq_in_i[i]+Delta_u  ##  Acción de control dentro de los límites establecidos, esta se actualiza normalmente sumandole el cambio.  

    if Nivel_i[i]>=Diam*0.75:  ##  Corrige el nivel cuando se alcanza su máximo permitido.
        Nivel_i[i]=Diam*0.75
        print("Nivel e líquido máximo de permitido")  
        continue
        
    elif Nivel_i[i]<0:  ##  Corrige el nivel cuando se alcanza su mínimo permitido.
        Nivel_i[i]=0
        print("Nivel de líquido mínimo  permitido")  
        continue
    
    else:
        Nivel_i[i]=Nivel  ##  El nivel esta dentro de los límites pertimitos.
        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.6.2-  ACCIÓN DE CONTROL PARA LA PRESIÓN DE VAPOR. 
          
    if (Q_i[i]+Delta_uP)>Q_Max: 
        Q_i[i+1]=Q_Max  ## Corrige la acción de control para que no pase de su límite máximo.
        
    elif (Q_i[i]+Delta_uP)<0:
        Q_i[i+1]=0  ##  Corrige la acción de control para que no pase de su límite mínimo.
        
    else:
        Q_i[i+1]=Q_i[i]+Delta_uP  ##  Acción de control dentro de los límites establecidos, se actualiza normalmente sumandole el cambio.   

    if P_Vap_i[i]<=0.04:
        P_Vap_i[i]=0.04  ##  Corrige la presión de vapor cuando se alcanza su mínimo permitido.
        print("Presión de vapor mínima permitida")
        continue
        
    elif P_Vap_i[i]>150:  ##  Corrige la presión de vapor cuando se alcanza su máximo permitido.
        P_Vap_i[i]=150
        print("Presión de vapor máxima permitida")
        continue
    
    else:
        P_Vap_i[i]=P_Vapor  ##  La presión de vapor esta dentro de los límites pertimitos.
        
FlujoM_Liq_in_i=FlujoM_Liq_in_i[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
Q_i=Q_i[:i+1]   ##  Ajusta el vector a las dimensiones requeridas. 
  
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 10-  GRAFICAS Y RESULTADOS.        

pt.figure("MODELO CALDERA (CONTROL PID)", [12,7])  ##  Crea figura y le otorga dimensiones.
 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.1-  PRESIÓN DE VAPOR VS TIEMPO.

pt.subplot(3,3,1)  ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot) #|4 5 6|#.
pt.tight_layout(pad=5, w_pad=7, h_pad=4)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(tiempo, P_Vap_i, "b", linewidth=2)  ##  Grafica los datos requeridos (x , y , Color de la linea , Grosor de la linea).
pt.plot(tiempo, SP_P_i, "y--", label="Set point", linewidth=1.5)  ##  Grafica el Set Point.
pt.grid(True)  ##  Agrega la cuadrilla a la grafica.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel("Presión \nde vapor [Bar]")  ##  Agrega título al eje y.
pt.title("Variables de proceso")
pt.ylim(min(P_Vap_i)-0.1,max(P_Vap_i)+0.1)  ##  Ajusta los limites en el eje (y) para mejor visualización.
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
pt.legend(loc="best")  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.2-  GRAFICA LOS DATOS, FLUJO MÁSICO DE VAPOR QUE SALE DE LA CALDERA VS TIEMPO.

pt.subplot(3, 3, 2)
pt.plot(tiempo, FlujoM_Vap_out_i, "g", linewidth=2)
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Flujo de vapor\n salida [kg/s]") 
pt.title("Perturbaciónes")
pt.ylim(min(FlujoM_Vap_out_i)-0.1,max(FlujoM_Vap_out_i)+0.1) 
pt.tick_params(labelsize=8) 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.3-  FLUJO DE ENERGÍA VS TIEMPO.
pt.subplot(3, 3, 3)
pt.plot(tiempo, Q_i, "b", linewidth=2) 
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Flujo de\n calor [kw]") 
pt.ylim(min(Q_i)-100,max(Q_i)+100)
pt.title("Variables de control") 
pt.tick_params(labelsize=8)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.4-  NIVEL VS TIEMPO.

pt.subplot(3, 3, 4) 
pt.plot(tiempo, Nivel_i, "k", linewidth=2)
pt.plot(tiempo, SP_N_i, "y--", label="Set point", linewidth=1.5)
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Nivel de \nlíquido [m] ")
pt.ylim(min(Nivel_i)-0.0001,max(Nivel_i)+0.0001)
pt.tick_params(labelsize=8)  
pt.legend(loc="best") 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.5-  FLUJO MÁSICO VAPOR INICIAL O LÍQUIDO QUE SE EVAPORA VS TIEMPO.

pt.subplot(3, 3, 5)
pt.plot(tiempo, FlujoM_LiqEvap_i, "y", linewidth=2) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Flujo de líquido\nevaporado [kg/s]") 
pt.ylim(min(FlujoM_LiqEvap_i)-0.1,max(FlujoM_LiqEvap_i)+0.1)  
pt.tick_params(labelsize=8) 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.6-  FLUJO MÁSICO DE LÍQUIDO VS TIEMPO.

pt.subplot(3, 3, 6)
pt.plot(tiempo, FlujoM_Liq_in_i, "c", linewidth=2) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Flujo de\n líquido [kg/s]") 
pt.ylim(min(FlujoM_Liq_in_i)-0.1,max(FlujoM_Liq_in_i)+0.1)  
pt.tick_params(labelsize=8) 
 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.7-  ERROR EN LA PRESIÓN DE VAPOR VS TIEMPO.
pt.subplot(3, 3, 7)
pt.plot(tiempo, Error_P_i, "r", linewidth=2) 
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Error Presión [Bar]") 
pt.title("Error de Presión")
pt.ylim(min(Error_P_i)-0.1,max(Error_P_i)+0.1) 
pt.tick_params(labelsize=8) 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.8-  ERROR EN EL NIVEL DEL LÍQUIDO VS TIEMPO.

pt.subplot(3, 3, 8)
pt.plot(tiempo, Error_N_i, "g", linewidth=2) 
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Error Nivel [m]") 
pt.title("Error de Nivel")
pt.ylim(min(Error_N_i)-0.0001,max(Error_N_i)+0.0001) 
pt.tick_params(labelsize=8) 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.9-  VOLUMEN TANTO DEL LÍQUIDO COMO DEL VAPOR DENTRO DE LA CALDERA VS TIEMPO.

pt.subplot(3, 3, 9)
pt.plot(tiempo, Vol_Vap_i, "r", linewidth=2,label="Vapor")
pt.plot(tiempo, Vol_Liq_i, "b", linewidth=2,label="Líquido")  
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Volumen \ncaldera [m\N{superscript three}]")
pt.title("Volumen Liq-Vap en la caldera")
pt.tick_params(labelsize=8) 
pt.legend(loc="best")

#________________________________________________________________________________________________________________________________________________________________________
