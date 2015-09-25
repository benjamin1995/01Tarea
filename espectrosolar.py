import numpy as np
import matplotlib.pyplot as plt
from numpy import *                # Esencial para trabajar con arreglos
from scipy import integrate
from astropy import *
from astropy.constants import c,au,h,k_B,sigma_sb


#PARTE1

data = np.loadtxt("sun_AM0.dat")
londa=data[:,0] *(10**-3)          # De nm a micron
flujo=data[:,1] *(10**3)           # De W/m^2 /micron a erg/seg*cm2 /micron
plt.plot(londa,flujo,'r')
plt.xlabel('Longitud de onda (micron) ')
plt.ylabel('Flujo solar (erg/seg*cm2)/micron')
plt.title('Flujo solar v/s Longitud de onda')
plt.grid(True)
plt.savefig('flujolongitud.png')
plt.show()                    #Favor,al hacer zoom se verifica que hay una longitud de onda que maximiza el flujo solar


#PARTE2
#Integral de flujo sobre londa con Metodo del trapecio version simple 1

i=0
integral2=0
while i<=1695:        # Recorre todos los datos
  fa=data[i,1]        # flujo i
  fb=data[i+1,1]      # flujo i+1
  a=data[i,0]         # longitud de onda i
  b=data[i+1,0]       # longitud de onda i+1
  deltalonda=(b-a)
  integral2+=(deltalonda/2)*(fa+fb)  # Sumando areas
  i+=1

print "El valor numerico para la integral2 con el metodo del trapecio version 1 es de : ",integral2 ,"[W/m2]"

# luminosidad solar= 4*pi*radiosolar^2 * sigma * T^4
# luminosidad solar= 1366*S=1366*4*pi*ua^2
lumsolar= 4*np.pi*au**2 *integral2
print "El valor numerico para la Luminosidad solar es de :",lumsolar,"Watt/m2"

#PARTE3

#Integral3 (integral de fnormalizada entre 0 y pi/2)
#Definicion integral de la funcion de planck (factor constante se utilizara mas adelante)

n=100            # puede ir variando
a=0.04598        #Como en 0 y pi medio la funcion se comporta singularmente, acotamos un poco el intervalo.
b=np.pi/2 -a     #Dado que la funcion converge en b si nos acercamos por la izquiera
y=np.linspace(a,b,n)                                                #Arreglo desde a hasta b
fnormalizada= (np.sin(y)**3/np.cos(y)**5)*(1/(np.e**(np.tan(y))-1)) #La integral de f.Planck con cambio de variable arcotangente
                                                                    #Diverge cuando nos acercamos a b por la derecha y converge para a

#METODO DEL TRAPECIO VERSION 2
i=1
integral3=0
while i<=(n-4):
  deltay=(b-a)/n
  integral3+=(deltay/2)*(fnormalizada[1]+fnormalizada[n-2]+2*fnormalizada[i+1])
  i+=1

#Omito tanto el primer como el ultimo termino porque se indefine fnormalizada
print "El valor numerico para la integral3 con el metodo del trapecio version 2 es de : ", integral3


T=5778  #Temperatura superficial efectiva solar en kelvin
Ifplanck= ((2*np.pi*h)/c**2)*((k_B*T/h)**4)*integral3
print "El valor numerico para la radiacion de un cuerpo negro con la temperatura efectiva solar es de :",Ifplanck,"K4"


#RADIO EFECTIVO
radioeff=sqrt((lumsolar/Ifplanck)/(4*np.pi))
print "El radio efectivo es:",radioeff , "metros" ,"lo que equivale a aproximadamente", 6.68458333*10**-12 *radioeff,"U.A(1 radio solar)"


#PARTE4

#TEORICO PARA PREGUNTA 3

#INTEGRATE.QUAD

x3 = lambda x: x**3/(np.exp(x)-1)    # funcion x3 que puedo evaluar en los puntos que este definida la funcion
quad3=integrate.quad(x3,0, np.inf)
print "El valor teorico para la integral3 con el metodo quad es de :",quad3[0],"con un error de", quad3[1]


#INTEGRATE.TRAPZ

ynew=np.linspace(a,b,100)  # ni 0 ni pi/2 ya que se indefine integral
trapz3=integrate.trapz((np.sin(ynew)**3/np.cos(ynew)**5)*(1/(np.e**(np.tan(ynew))-1)),ynew)
print "El valor teorico para la integral3 con el metodo del trapecio es : ",trapz3


#TEORICO PARA PREGUNTA 2

#ES UN ARREGLO AL QUE SE LE APLICA METODO DE TRAPECIO

trapz2=np.trapz(flujo,londa)
print "El valor teorico para la integral2 con el metodo del trapecio es : ", trapz2 , "[W/m2]"


#TIEMPO PARA PREGUNTA 3
 #TIEMPO QUAD3
     #Puse en ipython %timeit quad3 y dio 10000000 loops,best of 3: 37.2 nanoseg por loop
 #TIEMPO TRAPZ3
     #Puse en ipython %timeit trapz3 y dio 10000000 loops,best of 3: 37.6 nanoseg por loop
 #TIEMPO NUMERICO3
     #Puse en ipython %timeit integral3 y dio 10000000 loops,best of 3: 38.9 nanoseg por loop


#TIEMPO PARA PREGUNTA 2
  #TIEMPO TRAPZ2
    #Puse en ipython %timeit trapz2 y dio 10000000 loops,best of 3: 37.3 nanoseg por loop
  #TIEMPO NUMERICO2
    #Puse en ipython %timeit integral2 y dio 10000000 loops,best of 3: 37.9 nanoseg por loop
