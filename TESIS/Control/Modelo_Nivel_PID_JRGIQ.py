# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 23:14:27 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.


#________________________________________________________________________________________________________________________________________________________________________


*** C O N T R O L  D E  N I V E L  P A R A  U N  T A N Q U E  (P,PI,PD,PID) .***

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO***


- El tanque de nivel estudiado a continuación, es un recipiente en el cuál ingresa y sale líquido.

- La geometria utilizada para el tanque es cilíndrica vertical.

- La variable de proceso es el nivel del tanque. Este se mantiene regulando el flujo de salida del líquido.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***


-Fluido de proceso = Agua líquida a T=25 [°C].

-Variable de proceso (y)= Nivel de agua en el tanque [m]."Nivel_i".

-variable de control = Flujo volumétrico de salida [m^3/s]."FlujoV_out_i".

-Variable de perturbación = Flujo volumétrico de entrada [m^3/s]."FlujoV_Liq_in_i".

-Acción de control (u) = Apertura válvula de salida [%]."Val_Ctrl_i".

-Accion de perturbación (d) = Apertura de la válvula de entrada [%]."Val_Pert_i"


**PARA MAS DETALLER VER MANUAL DE USUARIO**
#________________________________________________________________________________________________________________________________________________________________________
"""
#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np  ##  Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Area=1  ##  Área transversal del tanque [m^2].
Area_Orf=1.87536e-4  ##  Área del orificio de descarga del tanque [m^2].
SP=1.63165  ##  Asigna un valor al set point.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARÁMETROS DEL LAS SUSTANCIAS Y CONSTANTES.

g=9.8  ##  Constante gravitacional [m/s^2].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES EN ESTADO ESTACIONARIO.(ini=Inicial, in=Entrada).

Nivel_ini=SP  ##  Nivel inicial del tanque [m].
Val_Pert_ini=50  ##  Apertura inicial de la válvula de entrada (perturbación) [%].
Val_Ctrl_ini=50  ##  Apertura inicial de la válvula de salida (acción de control) [%].
FlujoV_Liq_in=0.001060537  ##  Flujo volumétrico de entrada.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN,PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final =300#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
Paso =1#float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ##  Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
Iterac = N_Part  ##  Número de iteraciones, ese puede decir que es igual al número de particiones.
t_Pert = int(0.3 * N_Part)  ##  Tiempo en que se aplica la perturbación 20% del tiempo final de simulación.
t_SP_2=int(0.7*N_Part)  ##  Tiempo donde cambia el Set Point (OPCIONAL).
n=1#int(input("Ingrese el tiempo de muestreo del controlador"))

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

FlujoV_Liq_in_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de entrada, (Variable de perturbación).
Nivel_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el nivel del tanque (y),(Variable proceso).
Val_Pert_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la apertura de la valvula entrada (d), (Acción de perturbación).
FlujoV_Liq_out_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de salida de la válvula (Variable control).
Error_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el error= Set Point- nivel.
SP_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point.
Area_Orf_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el coeficiente de la válvula.
Val_Ctrl_i = np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para la apertura de la válvula salida (u), (Acción de control).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 5-  SET POINT, PERTURBACIONES, Y VARIABLE DE CONTROL.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.1-  SET POINT.

SP_i[t_Inicial:]=SP  ##  Llena el vector del Set Point con el valor asignado como referencia.
#SP_i[t_SP_2:]=SP+0.01  ##  Cambia a un nuevo valor del Set Point, las posiciones del vector desde t_SP_2 en adelante.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.2-  PERTURBACIONES

Val_Pert_i[t_Inicial:]=Val_Pert_ini  ##  Inicia el valor de apertura para la válvula del líquido en estado estacionario [%].
Val_Pert_i[t_Pert:]=Val_Pert_ini+40  ##  Cambia el valor de apertura de la válvula.Tiempo en que inicia la perturbación t_Pert.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  PARÁMETROS DEL CONTROLADOR.

Kp=-1000  ##  Parámetro de ganancia proporcional [%].
ti=10  ##  Parámetro para el tiempo integral [s].
td=4  ##  Parámetro para el tiempo derivativo [s].
ITC=1  ##  Intervalo de tiempo de la acción de control.
e_1k=0  ##  Error en el instante de tiempo anterior (i-1).
e_2k=0  ##  Error en el instante de tiempo anterior (i-2).

#________________________________________________________________________________________________________________________________________________________________________
## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA DEL MODELO DEL PROCESO CON SU CONTROLADOR.##

Nivel=Nivel_ini  ##  Asigna el valor inicial al nivel del tanque para iniciar el método.

for i in range(0,Iterac,1):
    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 7-  ACTUALIZACIÓN DE LAS VARIABLES DE PROCESO, CONTROL Y ERROR.

    Nivel_i[i]= Nivel  ##  Guarda los datos para cada solución de la Ecuación diferencial del nivel.
    Error_i[i]=SP_i[i]-Nivel_i[i]  ##  Guarda el error en cada iteración.
    Val_Ctrl_i[t_Inicial]=Val_Ctrl_ini  ##  Guarda los datos para la variable de control.(Apertura de la válvula).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 8-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.1-  SOLUCIÓN DE ECUACIONES CONSTITUTIVAS.
    
    FlujoV_Liq_in_i[i]=(Val_Pert_i[i]/100)*FlujoV_Liq_in  ##  Soluciona y actualiza los valores para la ecuación de flujo volumétrico de entrada [m3/s].
    FlujoV_Liq_out_i[i]=(Val_Ctrl_i[i]/100)*Area_Orf*np.sqrt(2*g*Nivel_i[i])  ##  Soluciona y actualiza los valores para el flujo volumétrico de salida [m^3/s].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2-  SOLUCIÓN DE ECUACIONES DIFERENCIALES APLICANDO EL MÉTODO EULER.

    dLdt=(FlujoV_Liq_in_i[i]-FlujoV_Liq_out_i[i])/Area  ##  Soluciona la derivada de la ED para el nivel.
    Nivel=Nivel_ini+dLdt*Paso  ##  Soluciona la ecuación diferencial, agregandole al nivel inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    Nivel_ini=Nivel  ##  Actualiza el nivel para iniciar con el siguiente paso [m].
         
#________________________________________________________________________________________________________________________________________________________________________

#  MODULO 9-  SELECCIÓN Y ACCIONES DEL CONTROLADOR.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.1-  CONTROLADOR (P).

    if (td==0 and ti>=9999):  ##  Control (P), depende de los parametros td y ti.
        Delta_u = Kp*ITC*Error_i[i]  ##  Ecuación para el control (P).
 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.2-  CONTROLADOR (PI),(PD) Y (PID).

    else:                
         Delta_u= Kp*ITC*(Error_i[i]-e_1k) + Kp*(ITC/ti)*e_1k + Kp*(td/ITC)*(Error_i[i]-2*e_1k+e_2k)  ##  Ecuación para el control (PID)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.3-  ACTUALIZACION DEL ERROR.
         
    e_1k=Error_i[i]  ##  Actualiza el error anterior (i-1).
    e_2k=e_1k  ##  Actualiza el error anterior (i-2).
      
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.4-  ACCIÓNES DE CONTROL.    
    
    if (Val_Ctrl_i[i]+Delta_u)>100: 
        Val_Ctrl_i[i+1]=100  ##  Corrige la acción de control para que no pase de su límite máximo.
        
    elif (Val_Ctrl_i[i]+Delta_u)<0:        
        Val_Ctrl_i[i+1]=0  ##  Corrige la acción de control para que no pase de su límite mínimo.
    
    else:
         Val_Ctrl_i[i+1]=Val_Ctrl_i[i]+Delta_u  ##  Acción de control dentro de los límites establecidos, se actualiza normalmente sumandole el cambio.  

Val_Ctrl_i=Val_Ctrl_i[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
    
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 10-  GRAFICAS Y RESULTADOS.        

pt.figure("MODELO TANQUE (CONTROL PID)", [10,5])  ##  Crea figura y le otorga dimensiones.
 
#----------------------------------------------------------------------------------------------------------
#  MODULO 10.1-  GRAFICA LOS DATOS, NIVEL VS TIEMPO .
                                                                                                        #_______#
                                                                                                        #|1 2 3|#
pt.subplot(2, 3, 1)  ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot) #|4 5 6|#.
pt.tight_layout(pad=3, w_pad=5, h_pad=2)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.grid(True)  ##  Agrega la cuadrilla a la grafica.
pt.plot(tiempo, Nivel_i, "k", linewidth=3)  ##  Grafica los datos requeridos (x , y , Color de la linea , Grosor de la linea).
pt.plot(tiempo, SP_i, "y--", label="Set point", linewidth=2)  ##  Grafica el Set Point.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel("Nivel [m]")  ##  Agrega título al eje y.
pt.ylim(min(Nivel_i)-0.01,max(Nivel_i)+0.01)  ##  Ajusta los limites en el eje (y) para mejor visualización.
pt.title("Variable de proceso")
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
pt.legend(loc="best")  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.2- GRAFICA LOS DATOS, ACCIÓN DE PERTURBACIÓN (d) VS TIEMPO .

pt.subplot(2, 3, 2)
pt.grid(True)
pt.plot(tiempo, Val_Pert_i, "b", linewidth=3)
pt.xlabel("Tiempo [s]")
pt.ylabel("Apertura válvula \nEntrada [%]")
pt.title("Perturbación")
pt.tick_params(labelsize=8)
pt.ylim(0, 110)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.3- GRAFICA LOS DATOS, ACCIÓN DE CONTROL (u) VS TIEMPO.

pt.subplot(2, 3, 3)
pt.grid(True)
pt.plot(tiempo, Val_Ctrl_i, "g", linewidth=3)
pt.xlabel("Tiempo [s]")
pt.ylabel("Apertura válvula\nSalida [%]")
pt.title("Variable de control")
pt.tick_params(labelsize=8)
pt.ylim(0, 110)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 11.4- GRAFICA LOS DATOS, ERROR VS TIEMPO .

pt.subplot(2, 3, 4)
pt.grid(True)
pt.plot(tiempo, Error_i, "r", linewidth=3)
pt.xlabel("Tiempo [s]")
pt.ylabel("Error [m]")
pt.tick_params(labelsize=8)
pt.ylim(min(Error_i)- 0.001, max(Error_i)+ 0.001)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 11.5- GRAFICA LOS DATOS, FLUJO VÁLVULA ENTRADA (Variable de perturbación) VS TIEMPO.

pt.subplot(2,3,5)
pt.grid(True)
pt.plot(tiempo,FlujoV_Liq_in_i, 'c',linewidth=3)
pt.ylabel("Flujo entrada [m\N{superscript three}/s]")
pt.xlabel("Tiempo [s]")
pt.tick_params(labelsize=8)
pt.ylim(min(FlujoV_Liq_in_i)-0.001,max(FlujoV_Liq_in_i)+0.001)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.6- GRAFICA LOS DATOS, FLUJO VÁLVULA SALIDA (Variable Ctrl) VS TIEMPO.

pt.subplot(2, 3, 6)
pt.grid(True)
pt.plot(tiempo,FlujoV_Liq_out_i,'y',linewidth=3)
pt.ylabel("Flujo salida [m\N{superscript three}/s]")
pt.xlabel("Tiempo [s]")
pt.tick_params(labelsize=8)
pt.ylim(min(FlujoV_Liq_out_i)-0.001,max(FlujoV_Liq_out_i)+0.001)

#________________________________________________________________________________________________________________________________________________________________________
