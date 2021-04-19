import numpy as np
import matplotlib.pyplot as plt
def solve_linear_system(a, b, g, N):
    h = (b-a)/N
    tri_diag = np.eye(N-1,N+1)-2*np.eye(N-1,N+1,k=1)+np.eye(N-1,N+1,k=2)
    A = np.concatenate((np.concatenate((np.matrix(np.append(1,np.zeros(N))), tri_diag), axis=0), np.matrix(np.append(np.zeros(N), 1))), axis=0)
    b = np.append(0, np.append(-h*h*g*np.ones(N-1), 0))
    u = np.linalg.solve(A,b)

    return u
x = np.linspace(0,1,101)
y = solve_linear_system(0, 1, -10, 100)
plt.plot(x,y)
plt.show()


test
