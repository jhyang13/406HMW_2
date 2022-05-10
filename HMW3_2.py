from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math

#define axis
fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax2 = plt.axes(projection='3d')
ax3 = plt.axes(projection='3d')
ax4 = plt.axes(projection='3d')
#define variables and functions
#parameters
Qv = 413000
Qb = 280000
Qc = 280000
T_M = 3271
n = 4.2
D_0v = math.pow(10, -5) * 1.2
D_0b = math.pow(10, -14) * 5.7
D_0c = math.pow(10, -23)
Mu_0 = math.pow(10, 4) * 6.12
Mu_Temp_Dep = -0.42 #Tm/Mu_0 *dMu/dT
b = math.pow(10, -23) * 2.86
k = math.pow(10, -23) * 1.38
R = 8.314
Omega = math.pow(10, -29) * 1.8
delta = math.pow(10, -9) * 5
d = math.pow(10, -6) * 5
A = 7.5 * math.pow(10, 5)
As = math.pow(math.pow(3, 0.5), n+1) * A
Ac = A
A2 = As
#Temperature
T = [i for i in range(1, 3272)]
#Mu: shear modulus
Mu = []
for i in T:
    Mu.append(Mu_0 * (1 + (i - 300) / T_M * Mu_Temp_Dep))
#sigma_s
sigma_s = []
i, j = 0, 0
k = math.pow(10, -8)
while i <= math.pow(10, -3):
    j = j + 1
    sigma_s.append(k)
    k = k + 9 * math.pow(10, -7)
    if j >= 3271:
        break
#D_v
D_v = []
for i in T:
    D_v.append(D_0v * math.exp(-Qv/(R * i)))
#D_b
D_b = []
for i in T:
    D_b.append(D_0b * math.exp(-Qb/(R * i)))
#D_c
D_c = []
for i in T:
    D_c.append(D_0c * math.exp(-Qc/(R * i)) / Ac)
#D_eff
D_eff = []
for i in range(0, 3271):
    D_eff.append(D_v[i] * (1 + 10 * Ac / math.pow(b, 2) * math.pow(sigma_s[i] / Mu[i], 2) * D_c[i] / D_v[i]))
    print(D_eff[i])
#high temperature power-law creep
rate_highTcreep = []
for i in range(0, 3271):
    rate_highTcreep.append(A2 * D_eff[i] * Mu[i] * b / (k * T[i])*(math.pow(sigma_s[i] / Mu[i], n)))
#low temperature power-law creep
rate_lowTcreep = []
for i in range(0, 3271):
    rate_lowTcreep.append(A2 * D_eff[i] * Mu[i] * b / (k * T[i])*(math.pow(sigma_s[i], n+2)/(math.pow(Mu[i],n))))
#boundary diffusional flow
rate_highTdiff = []
for i in range(0, 3271):
    rate_highTdiff = 42 * sigma_s[i] * Omega * 3.14 * delta / (k * T[i] * math.pow(d, 3)) * D_b[i]
#lattice diffusional flow
rate_lowTdiff = []
for i in range(0, 3271):
    rate_lowTdiff = 42 * sigma_s[i] * Omega / (k * T[i] * math.pow(d, 2))*D_v[i]

#plot
surf1 = ax1.plot_surface(T, sigma_s, rate_highTcreep, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf2 = ax2.plot_surface(T, sigma_s, rate_lowTcreep, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf3 = ax3.plot_surface(T, sigma_s, rate_highTdiff, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf4 = ax4.plot_surface(T, sigma_s, rate_lowTdiff, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.show()
