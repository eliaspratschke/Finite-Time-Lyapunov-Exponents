import numpy as np
import matplotlib.pyplot as plt

ftle = np.genfromtxt('ftlelong.txt', delimiter=' ')


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(ftle, interpolation='nearest')
fig.colorbar(cax)

ax.set_xticklabels([''])
ax.set_yticklabels([''])

plt.savefig('ftle.png', dpi = 800)
plt.show()
