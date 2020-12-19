# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 17:36:11 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.
"""



#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.
import math as mt
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Longitud= 3  ##  Longitud de la torre [m].
Poros= 0.75  ##  Porosidad del empaque.
CorreHoldup= 0.0333  ##  Corrección de porosidad por retenido del líquido.
PorosOper= Poros-CorreHoldup  ##  Porosidad operativa, se resta la porosidad por retenido de líquido.
DiamEmp= 0.0472  ##  Diámtro del empaque utilizado.
aM= 37.4  ##  Área de transferencia de masa m^2/m^3.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.

MR= 11  ##  Masa molar refinado R [kg/kmol].
ME= 260  ##  Masa molar extracto E [kg/kmol].
ViscosidadR= 1e-5  ##  Viscosidad de la fase refinado [kg/m s].
ViscosidadE= 2e-3  ##  Viscosidad de la fase Extracto [kg/m s].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES.(Estado estacionario).(ini=Inicial, in=Entrada, out=salida).


yA_in=0.4  ##  Concentración inicial del soluto en el refinado respecto a la corriente como un todo .
xA_in=0.05  ##  Concentración inicial del soluto en el extracto respecto a la corriente como un todo. 
yA_out=0.19  ##  Concentración final del soluto en el refinado respecto a la corriente como un todo.


YA_in=((yA_in)/(1-yA_in))  ##  Fracción molar de entrada del soluto en el refinado respecto al solvente libre de soluto.
XA_in=((xA_in)/(1-xA_in))  ##  Fracción molar de entrada del soluto en el extracto respecto al solvente libre de soluto.
YA_out=((yA_out)/(1-yA_out))  ##  Fracción molar de salida del soluto en el extracto respecto al solvente libre de soluto.

FluxMolR= 0.05  ##  Flux molar de entrada de la fase refinado [kmolMix/m^2 s].
FluxMolE= 0.065  ##  Flux molar de entrada de la fase extracto [kmolMix/m^2 s].

FluxE=FluxMolE*(1-XA_in)  ##  Flux del solvente puro en la fase extracto. 
FluxR=FluxMolR*(1-YA_in)  ##  Flux del solvente puro en la fase refinado. 

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3- LONGITUD DE LA TORRE Y NÚMERO DE PARTICIONES EN EL ESPACIO DE ESTA.

L_Inicial=0  ##  Valor de la longitud en el punto en donde inicia la simulación.
L_Final=Longitud  ##  Valor de la longitud final de la torre.
PasoEsp=0.05  ##  Valor para cada tramo en que se va a analizar la torre o paso espacial.
NumPasosEsp=int((L_Final-L_Inicial)/(PasoEsp))  ##  Número de pasos espaciales.
Longitud=np.linspace(L_Inicial,L_Final,NumPasosEsp)  ##  Crea un vector de longitud para poder graficarlo.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# #  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

ConcR=np.zeros(len(Longitud)+1)  ##  Crea un vector de datos para la concentración del refinado.
ConcE=np.zeros(len(Longitud)+1)  ##  Crea un vector de datos para la concentración del extracto.
MolesR=np.zeros(len(Longitud)+1) ##  Crea un vector de datos para las moles de refinado.
MolesE=np.zeros(len(Longitud)+1)  ##  Crea un vector de datos para las moles de extracto.

ViscR=np.zeros(len(Longitud))  ##  Crea un vector de datos para la viscosidad.
KY=np.zeros(len(Longitud))  ##  Crea un vector de datos para el coeficiente de TdeMasa.
ConcFict=np.zeros(len(Longitud))  ##  Crea un vector de datos para la concentración ficticia.
NA=np.zeros(len(Longitud))  ##  Crea un vector de datos para el flux de TdeMasa.
NA_tramo=np.zeros(len(Longitud))  ##  Crea un vector de datos para el flux de TdeMasa en cada partición de longitud.
ConcDeseada=np.zeros(len(Longitud))  ##  Crea un vector de datos para la concentración deseada (fines graficos).

ConcE_Fig=np.zeros(len(Longitud))  ##  Crea un vector de datos para la concentración del extracto en la CDE.
CDE=np.zeros(len(Longitud))  ##  Crea un vector de datos para cada punto de la CDE.
ConcR_LdeO=np.zeros(len(Longitud))  ##  Crea un vector de datos para la concentración del refinado en la linea de operación.
ConcE_LdeO=np.zeros(len(Longitud))  ##  Crea un vector de datos para la concentración de extracto en la linea de operación.


# #________________________________________________________________________________________________________________________________________________________________________

# # AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA.##

ConcR[0]=YA_in  ##  Otorga un valor inicial a la concentración de refinado para resolver el método iterativo.
ConcE[0]=XA_in  ##  Otorga un valor inicial a la concentración de extracto para resolver el método iterativo.
MolesR[0]=FluxMolR  ##  Otorga un valor inicial a las moles de refinado para resolver el método iterativo.
MolesE[0]=FluxMolE  ##  Otorga un valor inicial a las moles de extracto para resolver el método iterativo.

for i in range(0,NumPasosEsp,1):   
# #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
# #  MODULO 5-  SOLUCIÓN DE ECUACIONES.
# #------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# #  MODULO 5.1-  SOLUCIÓN DE ECUACIONES QUE REPRESENTAN EL MODELO.
    
    ViscR[i]=ViscosidadR+(ConcR[i]/100)**2  ##  Recalcula la viscosidad del refinado.
    KY[i]=MolesR[i]*0.15*((DiamEmp*MolesR[i]*MR)/(ViscR[i]*(1-PorosOper)))**-0.36  ##  Actualiza el coeficiente englobante de TdeMasa [kmol/m^2 s].
    
    
    ConcFict[i]=(9.524*ConcE[i]**3) -(7.302*ConcE[i]**2)+(2.064*ConcE[i])  ##  Polinomio que representa la concentración ficticia.
    
    NA[i]=KY[i]*(ConcR[i]-ConcFict[i])  ##  Recalcula el flux de TdeMasa [kmol/m^2 s].
    NA_tramo[i]=NA[i]*aM*PasoEsp  ##  Recalcula el flux de TdeMasa en cada partición de longitud [kmol/m^2 s].
    
    
    MolesR[i+1]=MolesR[i]-NA_tramo[i]  ##  Recalcula las moles en el  refinado [Kmol?m^2 s].
    MolesE[i+1]=MolesE[i]+NA_tramo[i]  ##  Recalcula las moles en el extracto [kmol?m^2 s].
    
    ConcR[i+1]=ConcR[i]-((NA_tramo[i])/(FluxR))  ##  Recalcula la concentración del refinado.
    ConcE[i+1]=ConcE[i]+((NA_tramo[i])/(FluxE))  ##  Recalcula la concentración del extracto.
    
    mLdeO=-FluxE/FluxR  ##  Pendiente de la linea de operación.
    

# #------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# #  MODULO 5.2-  ECUACIONES PARA LA CDE.
    
    PartCDE=YA_in/NumPasosEsp  ##  Valor de cada partición para graficar la CDE.
    PartLdeO=YA_out/NumPasosEsp  ##  Valor de cada partición para graficar la linea de operación.    
    ConcE_Fig[i]=i*PartCDE  ##  Construye punto a punto la concentracion del extracto para graficarlo en la CDE.
    CDE[i]=9.524*ConcE_Fig[i]**3 -7.302*ConcE_Fig[i]**2+2.064*ConcE_Fig[i]  ##  Polinomio que representa la CDE.


    ConcE_LdeO[i]=ConcE[0]+i*PartLdeO  ##  Aumenta a la concentración del extracto el valor de cada partición para graficarlo.
    ConcR_LdeO[i]=ConcR[0]+mLdeO*(ConcE_LdeO[i]-ConcE[0])  ##  Cálcula y actualiza la concentración del refinado en la CDE en función de la concentración del extracto.

    XA_final=ConcE[i+1]  ##  Otorga el valor a concentración final de A en el extracto para graficar dicho punto. 
    YA_final=ConcR[0]+mLdeO*(XA_final-ConcE[0])  ##  Cálcula la concentración final de A en el refinado para graficar dicho punto.
                
    ConcDeseada[i]=YA_out  ##  Llena el vector de concentración deseada con los datos de concentración de salida que se desea para poder graficar la linea recta.
       
ConcR=ConcR[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
ConcE=ConcE[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
MolesR=MolesR[:i+1]  ##  Ajusta el vector a las dimensiones requeridas.
MolesE=MolesE[:i+1]  ##  Ajusta el vector a las dimensiones requeridas.


# #________________________________________________________________________________________________________________________________________________________________________
# #  MODULO 6-  GRAFICAS Y RESULTADOS.        
# #------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# #  MODULO 6.1-  GRAFICA CONCENTRACIoNES DEL EXTRACTO Y REFINADO VS LONGITUD.        

pt.figure("TdeMasa Cocorriente", [10,5])  ##  Crea figura y le otorga dimensiones.
pt.subplot(1,2,1)   ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot).
pt.tight_layout(pad=5, w_pad=5, h_pad=4)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(Longitud,ConcR,'k',linewidth=2,label='Concentración en el Refinado')  ##  Grafica el longitud vs Concentración del refinado.
pt.plot(Longitud,ConcE,'b',linewidth=2,label='Concentración en el Extracto')  ##  Grafica el longitud vs Concentración del extracto.
pt.plot(Longitud,ConcFict,'y',linewidth=2,label='Concentración ficticia')  ##  Grafica el longitud vs Concentración ficticia.
pt.plot(Longitud,ConcDeseada,'r',linewidth=2)  ##  Grafica el longitud vs Concentración deseada.
pt.xlabel("Longitud [m]")  ##  Agrega título al eje x.
pt.ylabel('Concentración respecto\n del solvente libre de soluto ')  ##  Agrega título al eje y.
pt.title('Concentración vs longitud')  ##  Agrega título en la parte superior.
pt.grid(True)  ##  Agrega la cuadrilla a la imagen
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 7})
leg.get_frame().set_alpha(0.5)  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).
pt.show()

# #------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# #  MODULO 6.2-  GRAFICA DE LA CDE.  
      
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
