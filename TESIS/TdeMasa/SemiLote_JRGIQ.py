# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:15:18 2020

@author: Jheison René Gutiérrez Gómez. JRGIQ.
jhrgutierrezgo@unal.edu.co.

#________________________________________________________________________________________________________________________________________________________________________


*** S I M U L A D O R  D I N A M I C O  D E  U N  P R O C E S O  D E  T R A N S F E R E N C I A  D E  M A S A  E N  S E M I  L O T E  ***

#________________________________________________________________________________________________________________________________________________________________________
*** DESCRIPCIÓN DEL PROCESO***


- Se simulara el proceso de TdeMasa en una torre con un lecho poroso.

- Proceso en semilote en donde una sola de las fases se encuentra en movimiento y la otra permanece estática.

- La fase que se encuentra en movimiento es la fase refinado que se supondra como un gas que atravieza el lecho poroso 'FluxR'.

- La fase estacionaria es el refinado que se supondra como el sólido del empaque de la torre 'molesE'.

- Debe tenerse en cuenta que este simulador depende tanto del tiempo como del espacio en cada tramo del lecho de la torre.

#________________________________________________________________________________________________________________________________________________________________________
***  DATOS DEL PROCESO  ***

-  Fase refinado R (GAS)
-  Fase extracto E (SÓLIDO)

#________________________________________________________________________________________________________________________________________________________________________
"""
#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Diam_Emp=0.05  ##  Díametro del empaque usado [m].
Long=2  ##  Longitud de la torre [m].
FluxMol_R=0.065  ##  Flux  molar del refinado [kmol/m{2 s}].
CargaMol_E=1  ##  Carga molar del extracto [kmol].
Area_F=0.9  ##  Área de flujo de la torre [m^2].
aM=  37.4  ##  Área de tranferencia de masa [m^2].
Porosidad=0.55  ##  Porosidad del lecho.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.

MMolar_R=11  ##  Masa molar refinado R [kg/kmol].
MMolar_E=35  ##  Masa molar extracto E [kg/kmol].
Viscosidad=1e-5  ##  Viscosidad del refinado[kg/m s]. 1e-5
DensidadM_Emp=5  ##  Densidad del empaque de la torre [kmol/m^3].


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES.(Estado estacionario).(ini=Inicial, in=Entrada, out=salida).

yA_in=0.3  ##  Concentración inicial del soluto en el refinado respecto a la corriente como un todo.
xA_in=0.05  ##  Concentración inicial del soluto en el extracto respecto a la corriente como un todo. 
yA_out=0.18  ##  Concentración final del soluto en el refinado respecto a la corriente como un todo.

YA_in=((yA_in)/(1-yA_in))  ##  Fracción molar de entrada del soluto en el refinado respecto al solvente libre de soluto.
XA_in=((xA_in)/(1-xA_in))  ##  Fracción molar de entrada del soluto en el extracto respecto al solvente libre de soluto.
YA_out=((yA_out)/(1-yA_out))  ##  Fracción molar de salida del soluto en el extracto respecto al solvente libre de soluto.

CargaE=CargaMol_E*(1-XA_in)  ##  Carga de solvente puro en la fase extracto [kmol]. 
FluxMolR=FluxMol_R*(1-YA_in)  ##  Flux del solvente puro en la fase refinado [Kmol/m^2 s]. 
ConcSat_E=0.28  ##  Concentración de saturación del empaque.


#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN, PASO Y NÚMERO DE PARTICIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 3.1- TIEMPO DE SIMULACIÓN Y NÚMERO DE PARTICIONES EN EL TIEMPO.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final = 500#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
PasoTemp = 1 ##float(input("ingrese el tamaño del PasoTemp para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
NumPasosTemp=int((t_Final-t_Inicial)/(PasoTemp))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, NumPasosTemp)  ##  Crea vector de tiempo para poder graficarlo [s].


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 3.2- LONGITUD DE LA TORRE Y NÚMERO DE PARTICIONES EN EL ESPACIO DE ESTA.

L_Inicial=0  ##  Valor de la longitud en el punto en donde inicia la simulación.
L_Final=Long  ##  Valor de la longitud final de la torre.
PasoEsp=0.05  ##  Valor para cada tramo en que se va a analizar la torre o paso espacial.

PasoEspVol=PasoEsp*Area_F  ##  Paso espacial para cada tramo del volumen.
NumPasosEsp=int((L_Final-L_Inicial)/(PasoEsp))  ##  Número de pasos espaciales.
Longitud=np.linspace(L_Inicial,L_Final,NumPasosEsp)  ##  Crea un vector de longitud para poder graficarlo.
CargaTramo=DensidadM_Emp*PasoEspVol  ##  Valor de la carga molar del sólido para cada tramo de la torre.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 3.3- PARTICIONES PARA CONSTRUIR LA CDE Y LA LINEA DE OPERACIÓN.

DiscretCDE=YA_in/NumPasosTemp  ##  Discretización de la fracción molar de soluto en el refinado a la entrada para graficar en la CDE.
DiscretLdeO=YA_out/NumPasosTemp  ##  Discretización de la fracción molar de soluto en el refinado a la salida para graficar la linea de operación.


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

ConcR=np.zeros((len(tiempo),len(Longitud)+1))  ##  Crea un vector de datos para la concentración del refinado.
ConcE=np.zeros((len(tiempo)+1,len(Longitud)))  ##  Crea un vector de datos para la concentración del extracto.
FluxR=np.zeros((len(tiempo)+1,len(Longitud)+1)) ##  Crea un vector de datos para las moles de refinado.
MolesE=np.zeros((len(tiempo)+1,len(Longitud)+1))  ##  Crea un vector de datos para las moles de extracto.

Visc=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para la viscosidad.
K_Liq=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para el coeficiente de TdeMasa.
KY=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para el coeficiente de TdeMasa englobante.
ConcFict=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para la concentración ficticia.
NA=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para el flux de TdeMasa.
NA_DeltaL=np.zeros((len(tiempo),len(Longitud)))  ##  Crea un vector de datos para el flux de TdeMasa en cada partición de tiempo.
ConcDeseada=np.zeros((len(tiempo),len(Longitud)))   ##  Crea un vector de datos para la concentración deseada (fines graficos).

ConcE_Fig=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración del extracto en la CDE.
CDE=np.zeros(len(tiempo))  ##  Crea un vector de datos para cada punto de la CDE.
ConcR_LdeO=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración del refinado en la linea de operación.
ConcE_LdeO=np.zeros(len(tiempo))  ##  Crea un vector de datos para la concentración de extracto en la linea de operación.


#________________________________________________________________________________________________________________________________________________________________________
## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA.##

FluxR[0][0]=FluxMol_R  ##  Otorga un valor inicial a las moles de refinado para resolver el método iterativo.
MolesE[0][0]=CargaMol_E/NumPasosEsp  ##  Otorga un valor inicial a las moles de extracto para resolver el método iterativo.

for m in range(0,NumPasosTemp,1): 
    for n in range(0,NumPasosEsp,1): 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
#  MODULO 5-  SOLUCIÓN DE ECUACIONES.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------  
#  MODULO 5.1-  CREACIÓN DE DATOS INICIALES PARA LAS MATRICES DE CONCENTRACIÓN.
        
        ConcR[m][n]=YA_in  ##  Otorga un valor inicial a la concentración de refinado para resolver el método iterativo.
        ConcE[m][n]=XA_in  ##  Otorga un valor inicial a la concentración de extracto para resolver el método iterativo.
        
        ConcDeseada[m][n]=YA_out  ##  Otorga a la matriz de concentración deseada valores de la concentración de salida.(Fines graficos).


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 5.2-  SOLUCIÓN DE ECUACIONES QUE REPRESENTAN EL MODELO.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
for i in range(0,NumPasosTemp,1): 
    for j in range(0,NumPasosEsp,1):

        Visc[i][j]=Viscosidad+(ConcR[i][j]/100)**2  ##  Cálcula y actualiza la viscosidad del refinado.
        K_Liq[i][j]=FluxR[i][j]*1.15*((Diam_Emp*FluxR[i][j]*MMolar_R)/(Visc[i][j]*(1-Porosidad)))**-0.36  ##  Cálcula y actualiza el coeficiente englobante de TdeMasa [kmol/m^2 s].

    
        if (ConcSat_E-ConcE[i][j])>0:  ##  Condiciona y pone límites al coeficiente de TdeMasa englobante.
            
            KY[i][j]=K_Liq[i][j]*(ConcSat_E-ConcE[i][j])**0.5
        else:
            KY[i][j]=K_Liq[i][j]*0
                
        if ConcFict[i][j]>0.3:  ##  Condiciona o pone límites para la concentración en la concentración ficticia.
            ConcFict[i][j]=0.3
                    
        elif ConcFict[i][j]<1e-5:
            ConcFict[i][j]=1e-5
                
        else:
                    
            ConcFict[i][j]=(9.524*ConcE[i][j]**3) -(7.302*ConcE[i][j]**2)+(2.064*ConcE[i][j])  ##  Polinomio que representa la concentración ficticia.           
              
            
        NA[i][j]=KY[i][j]*(ConcR[i][j]-ConcFict[i][j])  ##  Cálcula y actualiza el flux de TdeMa2sa [kmol/m^2 s].
        NA_DeltaL[i][j]=NA[i][j]*aM*PasoEsp  ##  Cálcula y actualiza flujo molar en cada partición de tiempo [kmol/ s].
            
            
        FluxR[i][j+1]=FluxR[i][j]-NA_DeltaL[i][j]  ##  Recalcula el flux molar en el refinado refinado [Kmol].
            
        ConcR[i][j+1]=ConcR[i][j]-((NA_DeltaL[i][j])/(FluxMolR))  ##  Recalcula la concentración en el refinado.
                
        if ConcR[i][j+1]<0:  ##  Condiciona y pone límites a la concentración en el refinado.
            ConcR[i][j+1]=0  ##  Ya que la concentración no puede ser negativa esta se queda en cero.
                    
        MolesE[i+1][j+1]=MolesE[i][j]+NA_DeltaL[i][j]  ##  Recalcula las moles en el extracto [kmol].
        ConcE[i+1][j]=ConcE[i][j]+((NA_DeltaL[i][j])/(CargaTramo))  ##  Recalcula la concentración en el extracto.
        
        FluxR[i+1][0]=FluxR[0][0]  ##  Otorga nuevamente el valor inicial del flux molar en el refinado a la primer partición en todo el tiempo.
        MolesE[i+1][0]=MolesE[0][0]  ##  Otorga nuevamente el valor inicial a las moles en el extracto a la primer partición en todo el tiempo.  


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  MODULO 5.3-  ECUACIONES PARA LA CDE.
    
        mLdeO=((-CargaE/t_Final)/(FluxMolR)) ##  Pendiente de la linea de operación.
        ConcE_Fig[i]=i*DiscretCDE  ##  Construye punto a punto la concentracion del extracto para graficarlo en la CDE.
        CDE[i]=9.524*ConcE_Fig[i]**3 -7.302*ConcE_Fig[i]**2+2.064*ConcE_Fig[i]  ##  Polinomio para graficar la CDE, este se gráfica punto a punto de acuerdo a la discretización realizada para la CDE.
            
            
        ConcE_LdeO[i]=ConcE[0][0]+i*DiscretLdeO  ##  Aumenta a la concentración del extracto el valor de la partición asignada para graficarlo.
        ConcR_LdeO[i]=ConcR[0][0]+mLdeO*(ConcE_LdeO[i]-ConcE[0][0])  ##  Cálcula y actualiza la concentración del refinado en la CDE en función de la concentración del extracto.
            
        XA_final=ConcE[i+1]  ##  Otorga el valor a concentración final de A en el extracto. 
        YA_final=ConcR[0][0]+mLdeO*(XA_final-ConcE[0][0])  ##  Cálcula la concentración final de A en el refinado.

                 
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  GRÁFICAS Y RESULTADOS.        
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  Otorga valores para las particiones en el tiempo donde se desea ver el comportamiento gráfico del modelo.

tiempo1=NumPasosTemp*0  
tiempo2=int(NumPasosTemp*0.1)
tiempo3=int(NumPasosTemp*0.2)
tiempo4=int(NumPasosTemp*0.3)
tiempo5=int(NumPasosTemp*0.4)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#  Otorga un valores para las particiones en el espacio donde se desea ver el comportamiento gráfico del modelo.

parti1=NumPasosEsp*0
parti2=int(NumPasosEsp*0.1)
parti3=int(NumPasosEsp*0.2)
parti4=int(NumPasosEsp*0.5)
parti5=int(NumPasosEsp*0.9)


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# MODULO 6.1-  GRAFICA LA DISTRIBUCIÓN DE CONCENTRACIONES PARA EL REFINADO EN CADA PARTICIÓN DE LA TORRE EN INSTANTES DE TIEMPO DETERMINADOS.

pt.figure("TdeMasa Semi Lote. COMPORTAMIENTO DE LA CONCENTRACIÓN EN LA TORRE EN EL TIEMPO Y EL ESPACIO", [10,5])  ##  Crea figura y le otorga dimensiones.
pt.subplot(2,2,1)
pt.tight_layout(pad=5, w_pad=4, h_pad=4)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(Longitud,ConcR[tiempo1,1:],'r',label='Tiempo 0 [s]',linewidth=2)  ##  Grafica concentración en el refinado en función de la longitud a un instante de tiempo determinado.
pt.plot(Longitud,ConcR[tiempo2,1:],'k',label='Tiempo 50 [s]',linewidth=2)
pt.plot(Longitud,ConcR[tiempo3,1:],'b',label='Tiempo 100 [s]',linewidth=2)
pt.plot(Longitud,ConcR[tiempo4,1:],'y',label='Tiempo 150 [s]',linewidth=2)
pt.plot(Longitud,ConcR[tiempo5,1:],'g',label='Tiempo 200 [s]',linewidth=2)
pt.plot(Longitud,ConcDeseada[0][0:],'c',label='Concentración deseada',linewidth=2)
pt.grid(True)  ##  Agrega la cuadrilla a las gráficas.
pt.xlabel("Longitud [m]")  ##  Agrega el título al eje x.
pt.ylabel("Concentración en el Refinado")  ##  Agrega título al eje y.
pt.title("REFINADO (GAS)")  ##  Agrega título a la grafica.
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 6})  ##  Otorga propiedades al cuadro de leyendas.
leg.get_frame().set_alpha(0.5)  ##  Otorga la propiedade de que la leyenda se vea transparente.
pt.show()


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# MODULO 6.2-  GRAFICA LA DISTRIBUCIÓN DE CONCENTRACIONES PARA EL EXTRACTO EN CADA PARTICIÓN DE LA TORRE EN INSTANTES DE TIEMPO DETERMINADOS.

pt.subplot(2,2,2)
pt.plot(Longitud,ConcE[tiempo1][0:],'r', label='Tiempo = 0 [s]',linewidth=2)
pt.plot(Longitud,ConcE[tiempo2][0:],'k',label='Tiempo 50 [s]',linewidth=2)
pt.plot(Longitud,ConcE[tiempo3][0:],'b',label='Tiempo 100 [s]',linewidth=2)
pt.plot(Longitud,ConcE[tiempo4][0:],'y',label='Tiempo 150 [s]',linewidth=2)
pt.plot(Longitud,ConcE[tiempo5][0:],'g',label='Tiempo 200 [s]',linewidth=2)
pt.xlabel("Longitud [m]")
pt.ylabel("Concentración en el Extracto")
pt.title("EXTRACTO (SÓLIDO)")
pt.grid(True)
pt.tick_params(labelsize=8)  
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 6})
leg.get_frame().set_alpha(0.5)  


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# MODULO 6.3-  GRAFICA LA DISTRIBUCIÓN DE CONCENTRACIONES PARA EL REFINADO EN CADA PARTICIÓN DE TIEMPO EN DIFERENTES TRAMOS DE LA TORRE.

pt.subplot(2,2,3)
pt.plot(tiempo,ConcR[0:,parti1],'r',label='Longitud 0 [m]',linewidth=2)
pt.plot(tiempo,ConcR[0:,parti2],'k',label='Longitud 0.2 [m]',linewidth=2)
pt.plot(tiempo,ConcR[0:,parti3],'b',label='Longitud 0.4 [m]',linewidth=2)
pt.plot(tiempo,ConcR[0:,parti4],'y',label='Longitud 1 [m]',linewidth=2)
pt.plot(tiempo,ConcR[0:,parti5],'g',label='Longitud 1.8 [m]',linewidth=2)
pt.plot(tiempo,ConcDeseada[0:,0],'c',label='Concentración deseada',linewidth=2)
pt.grid(True)
pt.xlabel("tiempo [s]")
pt.ylabel("Concentración en el Refinado")
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 6})
leg.get_frame().set_alpha(0.5)  
pt.show()


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# MODULO 6.4-  GRAFICA LA DISTRIBUCIÓN DE CONCENTRACIONES PARA EL EXTRACTO EN CADA PARTICIÓN DE TIEMPO EN DIFERENTES TRAMOS DE LA TORRE.

pt.subplot(2,2,4)
pt.plot(tiempo,ConcE[1:,parti1],'r',label='Longitud 0 [m]',linewidth=2)
pt.plot(tiempo,ConcE[1:,parti2],'k',label='Longitud 0,2 [m]',linewidth=2)
pt.plot(tiempo,ConcE[1:,parti3],'b',label='Longitud 0.4 [m]',linewidth=2)
pt.plot(tiempo,ConcE[1:,parti4],'y',label='Longitud 1 [m]',linewidth=2)
pt.plot(tiempo,ConcE[1:,parti5],'g',label='Longitud 1.8 [m]',linewidth=2)
pt.grid(True)
pt.xlabel("tiempo [s]")
pt.ylabel("Concentración en el extracto")
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True, prop={'size': 6})
leg.get_frame().set_alpha(0.5)  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).
pt.show()


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
# MODULO 6.5-  GRAFICA LA CURVA DE DISTRIBUCIÓN DE EQUILIBRIO CDE.

pt.figure("TdeMasa Semi Lote. CURVA DE DISTRIBUCIÓN DE EQUILIBRIO CDE", [10,5])  ##  Crea figura y le otorga dimensiones.

pt.plot(ConcE_Fig,CDE,'k',label='CDE')
pt.plot(ConcE_LdeO,ConcR_LdeO,'y',label='Linea de operación')
pt.plot(XA_in,YA_in,'ko',label="Punto inicial de operación")
pt.plot(XA_final,YA_final,'ro',label='Punto final de operación')
pt.plot(ConcE[1:,6],ConcR[0:,6],'b',label="Evolución del proceso")
pt.xlabel('Concentración extracto')
pt.ylabel('Concentración refinado')
pt.title('Curva de distribución de equilibrio (CDE)')
pt.grid(True)
leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True,prop={'size': 7})
leg.get_frame().set_alpha(0.5)  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).
pt.show()


#________________________________________________________________________________________________________________________________________________________________________
