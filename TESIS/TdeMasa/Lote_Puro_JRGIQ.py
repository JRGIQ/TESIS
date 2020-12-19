# -*- coding: utf-8 -*-
"""
Created on Tue May 26 12:53:35 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.

#________________________________________________________________________________________________________________________________________________________________________


*** S I M U L A D O R  D I N A M I C O  D E  U N  P R O C E S O  D E  T R A N S F E R E N C I A  D E  M A S A  P O R  L O T E S  ***

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO***


- Se simulara un proceso de TdeMasa en una torre a la cual se le adicionara una carga molar de líquido y sólido.

- Proceso en lote puro donde la dos fases se cargan en la torre.

- La fase refinado se supondra como un liquido 'MolesR'.

- La fase extracto se supondra como un sólido de geometria esferica sin porosidad 'molesE'.

- Debe tenerse en cuenta que este simulador solo depende del tiempo mas no del espacio porque se supone que las concentraciones iniciales son iguales en toda la torre.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***

-  Fase refinado R
-  Fase extracto E

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

Diam_Part=0.002  ##  Díametro de partícula [m].
Vol_Part= (4/3)*mt.pi*(Diam_Part/2)**3  ##  Volumen de una partícula suponiendo geometría esferíca [m^3].
Area_Part=4*mt.pi*(Diam_Part/2)**2  ##  Área de la particula.

CargaMol_R=1.5  ##  Carga molar del refinado [kmol].
CargaMol_E=1.2  ##  Carga molar del extracto [kmol].



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.

MMolar_R=22  ##  Masa molar refinado R [kg/kmol].
MMolar_E=35  ##  Masa molar extracto E [kg/kmol].
Viscosidad=0.6e-3  ##  Viscosidad del refinado[kg/m s]

Densidad_R=1000  ##  Densidad del refinado [kg/m^3].
Densidad_E=1850   ##  Densidad del extracto [kg/m^3].

DensidadMol_R=Densidad_R/MMolar_R  ## Densidad molar refinado [Kmol/m^3].
DensidadMol_E=Densidad_E/MMolar_E  ##  Densidad molar extracto [kmol/m^3].
VolR=CargaMol_R*DensidadMol_R  ##  Volumen del refinado en el equipo [m^3].
VolE=CargaMol_E*DensidadMol_E  ##  Volumen del extracto en el equipo [m^3].
Num_Part=VolE/Vol_Part  ##  Número de partículas.
Area_M=Num_Part*Area_Part  ##  Área de tranferencia de masa [m^2].
Rel_Fases=VolE/VolR  ##  Relación de volumen de fases.


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES.(Estado estacionario).(ini=Inicial, in=Entrada, out=salida).


yA_in=0.3  ##  Concentración inicial del soluto en el refinado respecto a la corriente como un todo .
xA_in=0.05  ##  Concentración inicial del soluto en el extracto respecto a la corriente como un todo. 
yA_out=0.19  ##  Concentración final del soluto en el refinado respecto a la corriente como un todo.


YA_in=((yA_in)/(1-yA_in))  ##  Fracción molar de entrada del soluto en el refinado respecto al solvente libre de soluto.
XA_in=((xA_in)/(1-xA_in))  ##  Fracción molar de entrada del soluto en el extracto respecto al solvente libre de soluto.
YA_out=((yA_out)/(1-yA_out))  ##  Fracción molar de salida del soluto en el extracto respecto al solvente libre de soluto.

CargaE=CargaMol_E*(1-XA_in)  ##  Flux del solvente puro en la fase extracto. Recordar que estos flux son hipoteticos, ya que no ingresa ni sale masa del equipo. 
CargaR=CargaMol_R*(1-YA_in)  ##  Flux del solvente puro en la fase refinado. 

mLdeO=-CargaE/CargaR  ##  Pendiente de la linea de operación.


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN, PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final = 6000#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
Paso = 1 ##float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el número de iteraciones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
DiscretCDE=YA_in/N_Part  ##  Discretización de la fracción molar de soluto en el refinado a la entrada para graficar la CDE.
DiscretLdeO=YA_out/N_Part  ##  Discretización de la fraccióm molar de soluto en el refinado a la salida para graficar la linea de operación.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

ConcR=np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para la concentración del refinado.
ConcE=np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para la concentración del extracto.
MolesR=np.zeros(len(tiempo)+1) ##  Crea un vector de datos para las moles de refinado.
MolesE=np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para las moles de extracto.

Visc=np.zeros(len(tiempo))  ##  Crea un vector de datos para la viscosidad.
KY=np.zeros(len(tiempo))  ##  Crea un vector de datos para el coeficiente de TdeMasa.
ConcFict=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración ficticia.
NA=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flux de TdeMasa.
NA_deltat=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flux de TdeMasa en cada partición de tiempo.
ConcDeseada=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración deseada (fines gráficos).


ConcE_Fig=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración del extracto en la CDE.
CDE=np.zeros(len(tiempo))  ##  Crea un vector de datos para cada punto de la CDE.
ConcR_LdeO=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración del refinado en la linea de operación.
ConcE_LdeO=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración de extracto en la linea de operación.


#________________________________________________________________________________________________________________________________________________________________________

# AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA.##

ConcR[0]=YA_in  ##  Otorga un valor inicial a la concentración de refinado para resolver el método iterativo.
ConcE[0]=XA_in  ##  Otorga un valor inicial a la concentración de extracto para resolver el método iterativo.
MolesR[0]=CargaMol_R  ##  Otorga un valor inicial a las moles de refinado para resolver el método iterativo.
MolesE[0]=CargaMol_E  ##  Otorga un valor inicial a las moles de extracto para resolver el método iterativo.

for i in range(0,N_Part,1):   
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
#  MODULO 5-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 5.1-  SOLUCIÓN DE ECUACIONES QUE REPRESENTAN EL MODELO.
    
    Visc[i]=Viscosidad+(ConcR[i]/3)**2  ##  Recalcula la viscosidad del refinado.
    KY[i]=MolesR[i]*3.8e-9*((Diam_Part*MolesR[i]*MMolar_R)/(Visc[i]*(1-Rel_Fases)))**-0.33  ##  Actualiza el coeficiente englobante de TdeMasa [kmol/m^2 s].
    
    
    ConcFict[i]=(9.524*ConcE[i]**3) -(7.302*ConcE[i]**2)+(2.064*ConcE[i])  ##  Polinomio que representa la concentración ficticia.
    
    NA[i]=KY[i]*(ConcR[i]-ConcFict[i])  ##  Recalcula el flux de TdeMasa [kmol/m^2 s].
    NA_deltat[i]=NA[i]*Area_M*Paso  ##  Recalcula las moles que se transfieren en cada partición de tiempo [kmol].
    
    
    MolesR[i+1]=MolesR[i]-NA_deltat[i]  ##  Recalcula las moles en el  refinado [Kmol].
    MolesE[i+1]=MolesE[i]+NA_deltat[i]  ##  Recalcula las moles en el extracto [kmol].
    
    ConcR[i+1]=ConcR[i]-((NA_deltat[i])/(CargaR))  ##  Recalcula la concentración del refinado.
    ConcE[i+1]=ConcE[i]+((NA_deltat[i])/(CargaE))  ##  Recalcula la concentración del extracto.
    
    
    

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 5.2-  ECUACIONES PARA LA CDE.
          
    ConcE_Fig[i]=i*DiscretCDE  ##  Construye punto a punto la concentracion del extracto para graficarlo en la CDE.
    CDE[i]=9.524*ConcE_Fig[i]**3 -7.302*ConcE_Fig[i]**2+2.064*ConcE_Fig[i]  ##  Polinomio para graficar la CDE, este se gráfica punto a punto de acuerdo a la discretización realizada para la CDE.


    ConcE_LdeO[i]=ConcE[0]+i*DiscretLdeO  ##  Aumenta a la concentración del extracto el valor de cada partición para graficarlo.
    ConcR_LdeO[i]=ConcR[0]+mLdeO*(ConcE_LdeO[i]-ConcE[0])  ##  Cálcula y actualiza la concentración del refinado en la CDE en función de la concentración del extracto.

    XA_final=ConcE[i+1]  ##  Otorga el valor a concentración final de A en el extracto para graficar dicho punto. 
    YA_final=ConcR[0]+mLdeO*(XA_final-ConcE[0])  ##  Cálcula la concentración final de A en el refinado para graficar dicho punto.
                
    ConcDeseada[i]=YA_out  ##  Llena el vector de concentración deseada con los datos de concentración de salida que se desea para poder graficar la linea recta.
       
ConcR=ConcR[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
ConcE=ConcE[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
MolesR=MolesR[:i+1]  ##  Ajusta el vector a las dimensiones requeridas.
MolesE=MolesE[:i+1]  ##  Ajusta el vector a las dimensiones requeridas.


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  GRAFICAS Y RESULTADOS.        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 6.1-  GRAFICA CONCENTRACIoNES DEL EXTRACTO Y REFINADO VS TIEMPO DE OPERACIÓN.        

pt.figure("TdeMasa Lote puro", [10,5])  ##  Crea figura y le otorga dimensiones.
pt.subplot(1,2,1)   ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot).
pt.tight_layout(pad=5, w_pad=5, h_pad=4)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(tiempo,ConcR,'k',linewidth=2,label='Concentración en el Refinado')  ##  Grafica el tiempo vs Concentración del refinado.
pt.plot(tiempo,ConcE,'b',linewidth=2,label='Concentración en el Extracto')  ##  Grafica el tiempo vs Concentración del extracto.
pt.plot(tiempo,ConcFict,'y',linewidth=2,label='Concentración ficticia')  ##  Grafica el tiempo vs Concentración ficticia.
pt.plot(tiempo,ConcDeseada,'r',linewidth=2)  ##  Grafica el tiempo vs Concentración deseada.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel('Concentración respecto\n del solvente libre de soluto ')  ##  Agrega título al eje y.
pt.title('Concentración vs tiempo')  ##  Agrega título en la parte superior.
pt.grid(True)  ##  Agrega la cuadrilla a la imagen
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 7})
leg.get_frame().set_alpha(0.5)  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).
pt.show()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 6.2-  GRAFICA DE LA CDE.  
      
pt.subplot(1,2,2)
pt.grid(True)
pt.plot(ConcE_Fig,CDE,'k',label='CDE')
pt.plot(ConcE_LdeO,ConcR_LdeO,'y',label='Linea de operación guia (imaginaria)')
pt.plot(XA_in,YA_in,'ko',label="Punto inicial de operación")
pt.plot(XA_final,YA_final,'ro',label='Punto final de operación')
pt.plot(ConcE,ConcR,'b')
pt.xlabel('Concentración extracto')
pt.ylabel('Concentración refinado')
pt.title('Curva de distribución de equilibrio (CDE)')
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 7})
leg.get_frame().set_alpha(0.5)  
pt.show()
#________________________________________________________________________________________________________________________________________________________________________
