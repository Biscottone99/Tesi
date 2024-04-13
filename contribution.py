import matplotlib.pyplot as plt
import numpy as np
import math

# Definisci il nome del file
filename = 'input_2.dat'

# Apri il file in modalità lettura
with open(filename, 'r') as file:
    number_str = file.readline().strip()
    rad = float(number_str)
    rad = math.degrees(rad)
data = np.loadtxt('S1_2.dat')

# Estrai colonne da dati
x = data[:, 0]  # Prima colonna
y1 = data[:, 1]  # Seconda colonna
y2 = data[:, 2]  # Terza colonna
#y3 = data[:, 3]  # Quarta colonna
#y4 = data[:, 4]  # Quinta colonna

# Tronca x a tre cifre decimali
x_troncato = np.round(x, 3)

# Larghezza delle barre
bar_width = 0.2  # Riduci la larghezza delle barre

# Posizioni delle barre sul grafico
index = np.arange(len(x))

# Creazione del grafico
plt.bar(index , y1, bar_width, color='b', label='Mono real')
plt.bar(index + bar_width, y2, bar_width, color='g', label='Mono cplx')
#plt.bar(index + 0.5 * bar_width, y3, bar_width, color='orange', label='Tot real')
#plt.bar(index + 1.5 * bar_width, y4, bar_width, color='red', label='Tot cplx')

# Aggiunta delle etichette
plt.xlabel('Energy (eV)')
plt.ylabel('Energy (eV)')
plt.title('SOC contribution S1 Hubbard-model angle {}'.format(rad))
plt.xticks(index, x_troncato, rotation=45, ha='right')  # Ruota le etichette e posiziona a destra
plt.legend()

# Mostra il grafico
plt.tight_layout()
plt.show()
