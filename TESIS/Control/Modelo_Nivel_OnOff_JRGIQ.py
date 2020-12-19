# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 12:03:47 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.


#________________________________________________________________________________________________________________________________________________________________________


*** C O N T R O L  D E  N I V E L  P A R A  U N  T A N Q U E  (O N - O F F) .***

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO***


- El tanque de nivel estudiado a continuación, es un recipiente en el cuál ingresa y sale líquido.

- La geometria utilizada para el tanque es cilíndrica vertical.

- La variable de proceso es el nivel del tanque. Este se mantiene regulando el flujo de salida del líquido.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***


-Fluido de proceso = Agua líquida a T=25 [°C].

-Variable de proceso (y)= Nivel de agua en el tanque [m]."Nivel_i".

-variable de control = Flujo volumétrico de salida [m^3/s]."FlujoV_Liq_out_i".

-Variable de perturbación = Flujo volumétrico de entrada [m^3/s]."FlujoV_Liq_in_i".

-Acción de control (u) = Apertura válvula de salida [%]."Val_Ctrl_i".

-Accion de perturbación (d) = Apertura de la válvula de entrada [%]."Val_Pert_i"


**PARA MAS DETALLES VER MANUAL DE USUARIO**
#________________________________________________________________________________________________________________________________________________________________________
"""
#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np  ##  Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt  ##  Libreria de python que permite realizar gráficas.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Area=1  ##  Área transversal del tanque [m^2].
Area_Orf=1.87536e-4 * 2 ##  Área del orificio de descarga del tanque [m^2].
SP=1.63165  ##  Asigna un valor al set point.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DEL LAS SUSTANCIAS Y CONSTANTES.

g = 9.8  ##  Constante gravitacional [m/s^2].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES EN ESTADO ESTACIONARIO.(ini=Inicial, in=Entrada).

Nivel_ini=SP  ##  Nivel inicial del tanque [m].
Val_Pert_ini=50  ##  Apertura inicial de la válvula de entrada (perturbación) [%].
Val_Ctrl_ini=50  ##  Apertura inicial de la válvula de salida (acción de control) [%].
FlujoV_Liq_in=0.001060537  ##  Flujo volumétrico de entrada.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN,PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final =300 #int(input("Ingrese el tiempo de simulación [s] = "))  ## Tiempo en que finaliza la simulación [s].
Paso =1#float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ##  Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
Iterac = N_Part  ##  Número de iteraciones, ese puede decir que es igual al número de particiones.
t_Pert = int(0.3 * N_Part)  ##  Tiempo en que se aplica la perturbación 20% del tiempo final de simulación.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

FlujoV_Liq_in_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de entrada, (Variable de perturbación).
Nivel_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el nivel del tanque (y),(Variable proceso).
Val_Pert_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la apertura de la valvula entrada (d), (Acción de perturbación).
FlujoV_Liq_out_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de salida de la válvula (Variable control).
Error_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el error= Set Point- nivel.
SP_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point.
Val_Ctrl_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la apertura de la válvula salida (u), (Acción de control).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 5-  SET POINT, PERTURBACIONES, Y VARIABLE DE CONTROL.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.1-  SET POINT.

SP_i[t_Inicial:]=SP  ##  Llena el vector del Set Point con el valor asignado como referencia.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.2-  PERTURBACIONES

Val_Pert_i[t_Inicial:]=Val_Pert_ini  ##  Inicia el valor de apertura para la válvula del líquido en estado estacionario [%].
Val_Pert_i[t_Pert:]=Val_Pert_ini+40  ## Cambia el valor de apertura de la válvula.Tiempo en que inicia la perturbación t_Pert.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  CONDICIONES PARA LA ZONA MUERTA

Band_ZM=0.01  ##  (+-) Límite o banda para la zona muerta [m]
Zona_Muerta =0#int(input("SIN zona muerta marque (0)--CON zona muerta marque(1) = "))

#________________________________________________________________________________________________________________________________________________________________________

## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA DEL MODELO DEL PROCESO CON SU CONTROLADOR.##

Nivel=Nivel_ini  ##  Asigna el valor inicial al nivel del tanque para iniciar el método.

for i in range(0, Iterac, 1):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 7-  ACTUALIZACIÓN DE LAS VARIABLES DE PROCESO, CONTROL Y ERROR.

    Nivel_i[i] = Nivel  ##  Guarda los datos para cada solución de la Ecuación diferencial del nivel.
    Val_Ctrl_i[i] = Val_Ctrl_ini  ##  Guarda los datos para cada solución de la variable de control.(Apertura de la válvula).
    Error_i[i] = (SP_i[i] - Nivel)  ##  Guarda el error en cada iteración.
   
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 8-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.1-  SOLUCIÓN DE ECUACIONES CONSTITUTIVAS.
    
    FlujoV_Liq_in_i[i]=(Val_Pert_i[i]/100)*FlujoV_Liq_in  ##  Soluciona y actualiza los valores para la ecuación de flujo volumétrico de entrada [m3/s].
    FlujoV_Liq_out_i[i]=(Val_Ctrl_i[i]/100)*Area_Orf*np.sqrt(2*g*Nivel)  ##  Soluciona y actualiza los valores para el flujo volumétrico de salida [m^3/s].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2-  SOLUCIÓN DE ECUACIONES DIFERENCIALES APLICANDO EL MÉTODO EULER.
    
    dLdt=(FlujoV_Liq_in_i[i]-FlujoV_Liq_out_i[i])/Area  ##  Soluciona la derivada de la ED para el nivel.
    Nivel=Nivel_ini+dLdt*Paso  ##  Soluciona la ecuación diferencial, agregandole al nivel inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    Nivel_ini=Nivel  ##  Actualiza el nivel para iniciar con el siguiente paso [m].


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 9-  SELECCIÓN Y ACCIONES DEL CONTROLADOR.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.1-  APLICA PARA EL CONTROLADOR ON_OFF SIN ZONA MUERTA

    if Zona_Muerta==0:
        
        if Error_i[i]<0:  ##  Acción de control de regulación (corrección de la perturbacion).
            
            Val_Ctrl_ini=100  ##  Corrige la acción de control, cuando la variable de proceso es mayor al Set Point.
            
        else:
            
            Val_Ctrl_ini=0  ##  Corrige la acción de control, cuando la variable de proceso es menor al Set Point.
            
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.2-   APLICA PARA EL CONTROLADOR CON ZONA MUERTA
            
    else:  ##  Acción de control de regulación(corrección de la perturbacion).
        
        if abs(Error_i[i])<Band_ZM:
            
            Val_Ctrl_ini=Val_Ctrl_i[i]  ##  Mantiene la acción de control en el mismo punto.
            
        else:  ##  Control On-Off convencional.
            
            if Error_i[i]<0:
                
                Val_Ctrl_ini=100  ##  Corrige la variable de proceso, cuando esta es mayor al Set Point
                
            else:
                
                Val_Ctrl_ini=0  ##  Corrige la variable de proceso, cuando esta es menor al Set Point
                
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 10-  GRAFICAS Y RESULTADOS.        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.1- ESCOGE EL TIPO DE GRAFICO SELECCIONADO.

if Zona_Muerta==0:
    pt.figure("MODELO TANQUE (CONTROL ON-OFF SIN ZONA MUERTA)", [10,5]) 

else:
    pt.figure("MODELO TANQUE (CONTROL ON-OFF CON ZONA MUERTA)", [10,5]) 

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.2- GRAFICA LOS DATOS, NIVEL (y) (Variable de proceso) VS TIEMPO .
                                                                                                        #_______#
                                                                                                        #|1 2 3|#
pt.subplot(2, 3, 1)  ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot) #|4 5 6|#.
pt.tight_layout(pad=4, w_pad=5, h_pad=2)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.grid(True)  ##  Agrega la cuadrilla a la grafica.
pt.plot(tiempo, Nivel_i, "k", linewidth=2)  ##  Grafica los datos requeridos (x , y , Color de la linea , Grosor de la linea).
pt.plot(tiempo, SP_i, "y--", label="Set point", linewidth=2)  ##  Grafica el Set Point.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel("Nivel [m]")  ##  Agrega título al eje y.
pt.title("Variable de proceso")
pt.ylim(min(Nivel_i)-0.001,max(Nivel_i)+0.001)  ##  Ajusta los limites en el eje (y) para mejor visualización.
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
pt.legend(loc="best")  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.3- GRAFICA LOS DATOS, ACCIÓN DE PERTURBACIÓN (d) VS TIEMPO .

pt.subplot(2, 3, 2)
pt.grid(True)
pt.plot(tiempo, Val_Pert_i, "b", linewidth=2)
pt.xlabel("Tiempo [s]")
pt.ylabel("Apertura válvula \nEntrada [%]")
pt.title("Perturbación")
pt.tick_params(labelsize=8)
pt.ylim(0, 110)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.4- GRAFICA LOS DATOS, ACCIÓN DE CONTROL (u) VS TIEMPO.

pt.subplot(2, 3, 3)
pt.grid(True)
pt.plot(tiempo, Val_Ctrl_i, "g", linewidth=2)
pt.xlabel("Tiempo [s]")
pt.ylabel("Apertura válvula\nSalida [%]")
pt.title("Variable de control")
pt.tick_params(labelsize=8)
pt.ylim(0, 110)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.5- GRAFICA LOS DATOS, ERROR VS TIEMPO .

pt.subplot(2, 3, 4)
pt.grid(True)
pt.plot(tiempo, Error_i, "r", linewidth=2)
pt.xlabel("Tiempo [s]")
pt.ylabel("Error [m]")
pt.tick_params(labelsize=8)
pt.ylim(min(Error_i)- 0.001, max(Error_i)+ 0.001)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.6- GRAFICA LOS DATOS, FLUJO VÁLVULA ENTRADA (Variable de perturbación) VS TIEMPO.

pt.subplot(2,3,5)
pt.grid(True)
pt.plot(tiempo,FlujoV_Liq_in_i,'c',linewidth=3)
pt.ylabel("Flujo entrada [m\N{superscript three}/s]")
pt.xlabel("Tiempo [s]")
pt.tick_params(labelsize=8)
pt.ylim(min(FlujoV_Liq_in_i)-0.001,max(FlujoV_Liq_in_i)+0.001)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.7- GRAFICA LOS DATOS, FLUJO VÁLVULA SALIDA (Variable Ctrl) VS TIEMPO.

pt.subplot(2, 3, 6)
pt.grid(True)
pt.plot(tiempo,FlujoV_Liq_out_i,'y',linewidth=2)
pt.ylabel("Flujo salida [m\N{superscript three}/s]")
pt.xlabel("Tiempo [s]")
pt.tick_params(labelsize=8)
pt.ylim(min(FlujoV_Liq_out_i)-0.001,max(FlujoV_Liq_out_i)+0.001)

#________________________________________________________________________________________________________________________________________________________________________
 
