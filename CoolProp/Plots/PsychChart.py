from CoolProp.HumidAirProp import HAProps
import matplotlib.pyplot as plt
import numpy as np

fig=plt.figure(figsize=(10,8))
ax=fig.add_axes((0.15,0.15,0.8,0.8))
ax.set_xlim(-30,55)
ax.set_ylim(0,0.03)

Tdb = np.linspace(-30,55,100)+273.15
w=0*Tdb
for R in np.array([5,10,15,20,30,40,50,60,70,80,90,100])/100.0:
## for R in np.array([10,20])/100.0:
    for i in range(len(Tdb)):
        w[i]=HAProps('W','T',float(Tdb[i]),'P',101.325,'R',float(R))
    ax.plot(Tdb-273.15,w,'r')
        

plt.draw()
plt.show()