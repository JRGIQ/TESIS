# -*- coding: utf-8 -*-
"""
Created on Tue May  5 21:42:57 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.


#________________________________________________________________________________________________________________________________________________________________________


*** C O N T R O L  D E  T E M P E R A T U R A  P A R A  U N  I N T E R C A M B I A D O R  D E  C A L O R  (O N - O F F). 

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO. ***


- El intercambiador de calor estudiado a continuación, es el típico IdeCalor de tubos y coraza. Este opera con un fluido de servicio, el cual es vapor proveniente de 
    una caldera. A  este vapor se le puede modificar la Presion de saturación, pero en realidad la variable manipulada es el flujo másico de vapor que ingresa al IdeCalor.

- El fluido de proceso (fluido a calentar) es agua, esta ingresa aprox/ a 25°. El flujo que ingresa al IdeCalor es la variable que perturba el sistema, aunque tambien lo 
    podria ser la temperatura de este mismo, pero esta perturbación no posee un gran impacto en el sistema, ya que tendria que variar bastante para ver la perturbacion
    actuar en el sistema, y esto en la indutria no ocurre.
    
- El objetivo de control es asegurarse de que el fluido de proceso, al salir del IdeCalor tenga una temperatura aprox/ 55°C.

- Tanto las acciónes de control como de perturbación son realizadas por válvulas ON-OFF.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***


-Fluido de proceso = Agua líquida a T=(+-) 25 [°C] (Fluido por la coraza).

-Fluido de servicio = Vapor de agua proveniente de una caldera f(P_vapor), (Fluido por la tuberia).

-Variable de proceso (y)= Temperatura de salida del agua [°C]. "T_out_i".

-variable de control = Flujo másico del vapor suministrado al intercambiador de calor. "FlujoM_Vap_i".

-Variable de perturbación = Flujo volumétrico de entrada de agua [m^3/s]. "FlujoM_Liq_in_i".

-Acción de control (u) = Apertura válvula On-Off para de vapor [%]."Val_Ctrl_i".

-Accion de perturbación (d) =  Apertura de la válvula On-Off para el líquido [%]."Valv_Pert_i".


**PARA MAS DETALLER VER MANUAL DE USUARIO**
#________________________________________________________________________________________________________________________________________________________________________
"""
#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np  ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt  # Libreria de python que permite realizar graficas.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Volumen = 0.08  ##  Volumen del IdeCalor [m^3].
SP = 55+273.15  ##  Set Point [K].
P_Vapor=7  ##  Presión de vapor de la caldera [bar].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DEL LAS SUSTANCIAS Y CONSTANTES.

Densidad=1000  ##  Densidad agua[kg/m^3].
Masa=Volumen*Densidad  ##  Masa de agua dentro de la coraza [kg].
R=8.314  ##  Constante de los gases ideales [kj/kmol K].
M=18.02  ##  Masa molar agua [Kg/Kmol].

A=8.712  ##  Constante A para el Cp para el líquido. (Tabla C3. Smith Van Ness).
B=1.25*1e-03  ##  Constante B para el Cp para el líquido. (Tabla C3. Smith Van Ness).
C=-0.18*1e-06  ##  Constante C para el Cp para el líquido. (Tabla C3. Smith Van Ness).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES EN ESTADO ESTACIONARIO.(ini=Inicial, in=Entrada, out=Salida)

Val_Pert_ini = 50  ##  Apertura incial de la válcua de entrada (d) [%].
T_in = 20+273.15  ##  Temperatura de entrada (agua) al tanque [k].
Val_Ctrl_ini = 50  ## Apertura inicial de la válvula de vapor [%].
T_out_ini=SP  ##  Temperatura inicial de salida del proceso.
FlujoM_Liq_in=1 ##  Flujo másico de líquido que ingresa al IdeCalor [kg/s].
FlujoM_Vap_in=0.06879425144  ##  Flujo másico de vapor que ingresa al IdeCalor [kg/s].
FlujoM_Vap_Max=2*FlujoM_Vap_in
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN,PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final = 300#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
Paso =1 ##float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el núm5ero de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
Iter = N_Part  ##  Número de iteraciones, se puede decir que es igual al número de particiones.
t_Pert = int(0.3* N_Part)  ##  Tiempo y valor en que se aplica la perturbación [%].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

T_out_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la temperatura de salida, variable de proceso (y).
Valv_Pert_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la (acción) de apertura de la valvula de entrada de líquido (d).
Error_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el error= set point- T_out.
Val_Ctrl_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la (acción) de control o apertura del válvula de vapor (u).
SP_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point.
FlujoM_Liq_in_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de líquido (Fluido de proceso). 
Q_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el calor energía suministrada.
FlujoM_Vap_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo másico de vapor (Variable de control).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 5-  SET POINT, PERTURBACIONES, Y VARIABLE DE CONTROL.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.1-  SET POINT.

SP_i[t_Inicial:] = SP  ##  Llena el vector del Set Point con el valor asignado como referencia.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.2-  PERTURBACIONES

Valv_Pert_i[t_Inicial:] = Val_Pert_ini ##  Llena el vector para la perturbación con el valor inicial asignado.
Valv_Pert_i[t_Pert:] = Val_Pert_ini + 30  ##  Cambia el valor de la perturbación al vector a un tiempo t_Pert(tiempo perturbación).

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  CONDICIONES PARA LA ZONA MUERTA.

Band_ZM=1 ##  (+-) Límite o banda para la zona muerta [K].
Zona_Muerta = 0#int(input("SIN zona muerta marque (0)--CON zona muerta marque(1) = "))

#________________________________________________________________________________________________________________________________________________________________________
## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA DEL MODELO DEL PROCESO CON SU CONTROLADOR.##

T_out=T_out_ini

for i in range(0, Iter, 1):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 7-  ACTUALIZACIÓN DE LAS VARIABLES DE PROCESO, CONTROL Y ERROR.

    T_out_i[i] = T_out  ##  Guarda los datos para cada solución de la Ecuación diferencial de Temperatura.
    Val_Ctrl_i[i] = Val_Ctrl_ini  ##  Guarda los datos para cada solución de la variable de control (Apertura de la válvula de vapor).
    Error_i[i] = (SP_i[i] - T_out)  ##  Guarda el error en cada iteración.
    
    
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 8-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.1-  SOLUCIÓN DE ECUACIONES CONSTITUTIVAS.
    
    FlujoM_Liq_in_i[i] = (Valv_Pert_i[i] / 100) * FlujoM_Liq_in ##  Soluciona y actualiza la ecuación para el flujo de proceso (perturbación).
    FlujoM_Vap_i[i]=(Val_Ctrl_i[i]/100)*FlujoM_Vap_Max  ##  Soluciona y actualiza el flujo de vapor (acción de control).
    
    Cp_L=(A+(B*T_out)+(C*T_out**2))*((R)/(M))  ##  Soluciona y actualiza la ecuación para el Cp del líquido [Kj/Kg].(Tabla C3. Smith Van Ness).
    h_Vap=-0.0006*P_Vapor**3+0.1689*P_Vapor**2-20.951*P_Vapor+2284.4  ##  Cálcula el calor latente del vapor en función de la Presión de saturación.
    Q_i[i]=FlujoM_Vap_i[i]*h_Vap  ##  Soluciona la ecuación de balance termico, esta es la energía que entrega el vapor de la caldera al IdeCalor (Calor latente).
       
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2-  SOLUCIÓN DE ECUACIONES DIFERENCIALES APLICANDO EL MÉTODO EULER.
    
    dTdt =((FlujoM_Liq_in_i[i]/Masa)*(T_in-T_out_i[i]))+((1)/(Masa*Cp_L))*Q_i[i]  ##  Soluciona la derivada de la ED para la temperatura.
    T_out= T_out_ini + dTdt * Paso  ##  Soluciona la ecuación diferencial, agregandole a la temperatura inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    T_out_ini = T_out  ##  Actualiza la temperatura para iniciar con la siguiente partición. 
    
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 9-  SELECCIÓN DEL CONTROLADOR On-Off.    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.1-  APLICA PARA EL CONTROLADOR ON_OFF SIN ZONA MUERTA

    if Zona_Muerta==0:
        
        if Error_i[i]>=0:  ##  Acción de control de regulación (corrección de la perturbacion).
            
            Val_Ctrl_ini=100  ##  Corrige la acción de control, cuando la variable de proceso es menor al Set Point.
            
        else:
            
            Val_Ctrl_ini=0  ##  Corrige la acción de control, cuando la variable de proceso es mayor al Set Point.
                       
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.2-  APLICA PARA EL CONTROLADOR ON_OFF CON ZONA MUERTA

    else:  ##  Acción de control de regulación(corrección de la perturbacion).
        
        if abs(Error_i[i])<Band_ZM:
            
            Val_Ctrl_ini=Val_Ctrl_i[i]  ##  Mantiene la acción de control en el mismo punto.
            
        else:  ##  Control On-Off convencional.
            
            if Error_i[i]<0:
                
                Val_Ctrl_ini=0  ##  Corrige la variable de proceso, cuando esta es mayor al Set Point
                
            else:
                
                Val_Ctrl_ini=100  ##  Corrige la variable de proceso, cuando esta es menor al Set Point
                               
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 10-  GRAFICAS Y RESULTADOS.                        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.1-  ESCOGE EL TIPO DE GRAFICO SELECCIONADO.
                
if Zona_Muerta==0:  ##  Condición para que aparezca el título de (SIN ZONA muerta) en el recuadro grafico.
    
    pt.figure("MODELO INTERCAMBIADOR DE CALOR (CONTROL ON-OFF SIN ZONA MUERTA)", [10,5]) 

else:  ##  Condición para que aparezca el título de (CON ZONA muerta) en el recuadro grafico.
    
    pt.figure("MODELO INTERCAMBIADOR DE CALOR (CONTROL ON-OFF CON ZONA MUERTA)", [10,5])  
           
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.2-  GRAFICA LOS DATOS, TEMPERATURA DE SALIDA (y) (Variable de proceso) VS TIEMPO .
                                                                                                          #_______#
                                                                                                          #|1 2 3|#
pt.subplot(2, 3, 1)  ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot) #|4 5 6|#.
pt.tight_layout(pad=4, w_pad=5, h_pad=2)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(tiempo, T_out_i-273.15, "k", linewidth=2)  ##  Grafica los datos requeridos (x , y , Color de la linea , Grosor de la linea).
pt.plot(tiempo, SP_i-273.15, "y--", label="Set point", linewidth=3)  ##  Grafica el Set Point.
pt.grid(True)  ##  Agrega la cuadrilla a la grafica.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel("Temperatura\nSalida [°C] ")  ##  Agrega título al eje y.
pt.ylim(min(T_out_i-273.15)-0.1,max(T_out_i-273.15)+0.1)  ##  Ajusta los limites en el eje (y) para mejor visualización.
pt.title("Variable de proceso")
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
pt.legend(loc="best",prop={'size': 7})  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.3-  GRAFICA LOS DATOS, ACCIÓN DE PERTURBACIÓN (d) VS TIEMPO .

pt.subplot(2, 3, 2) 
pt.plot(tiempo, Valv_Pert_i, "r", linewidth=3)
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Apertura válvula\nLíquido [%]")
pt.title("Perturbación") 
pt.ylim(min(Valv_Pert_i)-10,max(Valv_Pert_i)+10)  
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.4-  GRAFICA LOS DATOS, ACCIÓN DE CONTROL (u) VS TIEMPO .

pt.subplot(2, 3, 3)
pt.plot(tiempo, Val_Ctrl_i, "b", linewidth=1) 
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Apertura válvula\n Vapor [%]") 
pt.ylim(min(Val_Ctrl_i)-10,max(Val_Ctrl_i)+10) 
pt.title("Variable de control")
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.5-  GRAFICA LOS DATOS, POTENCIA O CALOR TRANSFERIDO VS TIEMPO.

pt.subplot(2, 3, 4)
pt.plot(tiempo, Error_i, "c", linewidth=2) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Error [°C]") 
pt.ylim(min(Error_i)-1,max(Error_i)+1)  
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.6-  GRAFICA LOS DATOS, FLUJO FLUIDO DE PROCESO (Variable de perturbación) VS TIEMPO.

pt.subplot(2, 3, 5)
pt.plot(tiempo, FlujoM_Liq_in_i, "y", linewidth=3)
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Flujo agua [kg/s]") 
pt.ylim(min(FlujoM_Liq_in_i)-0.1,max(FlujoM_Liq_in_i)+0.1) 
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.7-  GRAFICA LOS DATOS, FLUJO DE VAPOR (Variable Ctrl) VS TIEMPO.

pt.subplot(2, 3, 6)
pt.plot(tiempo, FlujoM_Vap_i, "g", linewidth=1) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Flujo de vapor [kg/s]") 
pt.ylim(min(FlujoM_Vap_i)-0.01,max(FlujoM_Vap_i)+0.01)  
pt.tick_params(labelsize=8) 

#__________________________________________________________________________________________________________________#________________________________________________________________________________________________________________________________________________________________________
