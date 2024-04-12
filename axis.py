import numpy as np
import matplotlib.pyplot as plt

# Lettura dei dati da file
data = np.loadtxt('mcrot-tot-real.dat')
x = data[:,0]
mxr = data[:,1]
myr=data[:,2]
mzr=data[:,3]

data = np.loadtxt('mcrot-tot-cplx.dat')
x1=data[:,0]
mxi=data[:,1]
myi=data[:,2]
mzi=data[:,3]

data= np.loadtxt('bcrot-tot-real.dat')
x2=data[:,0]
bxr=data[:,1]
byr=data[:,2]
bzr=data[:,3]

data= np.loadtxt('bcrot-tot-cplx.dat')
x3=data[:,0]
bxi=data[:,1]
byi=data[:,2]
bzi=data[:,3]

#PLOT1
plt.subplot(2, 3, 1)
plt.plot(x, mxr, marker='o', color='blue', label='Mono-x-real')
plt.plot(x, mxi, marker='o', color='green', label='Mono-x-cplx')
plt.xlim(0,360)

plt.title('SOC-mono along x axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()


#PLOT2
plt.subplot(2, 3, 4)
plt.plot(x, bxr, marker='o', color='blue', label='Bi-x-real')
plt.plot(x, bxi, marker='o', color='green', label='Bi-x-cplx')
plt.xlim(0,360)

plt.title('SOC-bi along x axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()


#PLOT3
plt.subplot(2, 3, 2)
plt.plot(x, myr, marker='o', color='blue', label='Mono-y-real')
plt.plot(x, myi, marker='o', color='green', label='Mono-y-cplx')
plt.xlim(0,360)

plt.title('SOC-mono along y axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()


#PLOT4
plt.subplot(2, 3, 5)
plt.plot(x, byr, marker='o', color='blue', label='Bi-y-real')
plt.plot(x, byi, marker='o', color='green', label='Bi-y-cplx')
plt.xlim(0,360)

plt.title('SOC-bi along y axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()

#PLOT5
plt.subplot(2, 3, 3)
plt.plot(x, mzr, marker='o', color='blue', label='Mono-z-real')
plt.plot(x, mzi, marker='o', color='green', label='Mono-z-cplx')
plt.xlim(0,360)

plt.title('SOC-mono along z axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()


#PLOT6
plt.subplot(2, 3, 6)
plt.plot(x, bzr, marker='o', color='blue', label='Bi-z-real')
plt.plot(x, bzi, marker='o', color='green', label='Bi-z-cplx')
plt.xlim(0,360)

plt.title('SOC-bi along z axis vs torsion angle for GS-T_CT D-A system')
plt.xticks(np.arange(0, 361, 45))
plt.xlabel('Angle (Degrees)')
plt.ylabel('Energy(eV)')
plt.legend()


# Mostra il grafico
plt.legend()
plt.tight_layout()  # Opzionale: migliora la disposizione dei subplot
plt.show()
