import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#fig = plt.figure(figsize=(11.,8.))      # figura (ancho, largo)
ax = plt.subplot()      # subfigura

# Datos
data = np.loadtxt('energias_planetas.dat')
tiempo = data[:,0] # cojo el tiempo de la primera columna
energia = data[:,1]

# configurar ejes
ax.set_ylabel('Energía', fontname='DejaVu Sans', fontsize='12')
ax.set_xlabel('Tiempo', fontname='DejaVu Sans', fontsize='12')
ax.set_ylim(-0.0002,0.0002)

#Cambiar ticks
#for label in ax.get_xticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.xticks(fontsize='15')
#for label in ax.get_yticklabels():
    #label.set_fontproperties('Times New Roman')
#plt.yticks(fontsize='15')

# Creación de la gráfica
ax.plot(tiempo, energia, linestyle='-', marker='', markersize=4, color='#B4045F')  #marker=puntos

# Guardar la gráfica
plt.savefig('grafica_energia.png',dpi=300, bbox_inches = "tight")

# Mostrarla en pantalla
#plt.show()