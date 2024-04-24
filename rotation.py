import numpy as np

def leggi_coordinate(filename):
    return np.loadtxt(filename, usecols=(1, 3, 4, 5))

def traslazione(coord, traslazione):
    return coord + traslazione

def rotazione_xy(coord, angolo):
    angolo_rad = np.deg2rad(angolo)
    R = np.array([[np.cos(angolo_rad), -np.sin(angolo_rad), 0],
                  [np.sin(angolo_rad), np.cos(angolo_rad), 0],
                  [0, 0, 1]])
    return np.dot(R, coord.T).T

def rotazione_yz(coord, angolo):
    angolo_rad = np.deg2rad(angolo)
    R = np.array([[1, 0, 0],
                  [0, np.cos(angolo_rad), -np.sin(angolo_rad)],
                  [0, np.sin(angolo_rad), np.cos(angolo_rad)]])
    return np.dot(R, coord.T).T

def assegna_lettera(numero):
    if numero == 1:
        return 'H'
    elif numero == 6:
        return 'C'
    elif numero == 7:
        return 'N'
    elif numero == 8:
        return 'O'
    else:
        return 'X'  # Numero non valido, assegna 'X' come valore di default

def scrivi_coordinate(filename, coord, lettere):
    with open(filename, 'w') as file:
        for i in range(len(coord)):
            file.write(f"{lettere[i]} {coord[i, 0]:.5f} {coord[i, 1]:.5f} {coord[i, 2]:.5f}\n")

if __name__ == "__main__":
    # Lettura delle coordinate e numeri della seconda colonna da file
    coordinate_e_numeri = leggi_coordinate("coord.dat")
    numeri_colonna2 = coordinate_e_numeri[:, 0]
    coordinate = coordinate_e_numeri[:, 1:]
    
    # Traslazione delle coordinate
    traslazione_vettore = np.array([0.437975, 0.683213, -1.749063])  # Modifica con il vettore di traslazione desiderato
    coordinate_trasl = traslazione(coordinate, traslazione_vettore)
    
    # Rotazione nel piano xy di 342 gradi
    coordinate_rot_xy = rotazione_xy(coordinate_trasl, 342.24765)
    
    # Rotazione nel piano yz di 331.17787 gradi
    coordinate_rot_yz = rotazione_yz(coordinate_rot_xy, 331.17787)
    
   
    # Assegnazione delle lettere in base ai numeri della seconda colonna
    lettere = [assegna_lettera(numero) for numero in numeri_colonna2]
    
    # Scrivi le coordinate e le lettere nel file di output
    scrivi_coordinate("coordinate_trasformate.dat", coordinate_rot_yz, lettere)
