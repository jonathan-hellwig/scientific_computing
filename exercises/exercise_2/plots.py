import matplotlib
import matplotlib.pyplot as plt
import numpy as np

t = 1e5
k = np.arange(1,6)
p = np.power(10, k) 
S = t / (t/p + p / 10) 

fig, ax = plt.subplots()
ax.plot(k, S)

ax.set(xlabel='Number of rings', ylabel='Effective speedup',
       title='Effective Speedup')
ax.grid()

fig.savefig("speedup.png")
plt.show()

