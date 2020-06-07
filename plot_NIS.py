import numpy as np
import matplotlib.pyplot as plt



X, Y = [], []
for line in open('NIS_Radar.txt', 'r'):
  values = [float(s) for s in line.split()]
  X.append(values[0])
  
print(np.mean(X))
print(np.var(X))
plt.plot(X)
plt.show()