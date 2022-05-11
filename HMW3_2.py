import matplotlib.pyplot as plt
import math
import numpy as np

#define parameters
Qv = 413000
Qb = 280000
Qc = 280000
T_M = 3271
n = 4.2
D_0v = math.pow(10, -5) * 1.2
D_0b = math.pow(10, -14) * 5.7
D_0c = math.pow(10, -23)
Mu_0 = math.pow(10, 4) * 6.12
Mu_Temp_Dep = -0.42
b = math.pow(10, -23) * 2.86
R = 8.314
Omega = math.pow(10, -29) * 1.8
delta = math.pow(10, -9) * 5
d = math.pow(10, -6) * 5
A = 7.5 * math.pow(10, 5)
As = math.pow(math.pow(3, 0.5), n+1) * A
Ac = A
A2 = As

#define functions
#Temperature
T = np.linspace(100, 3272, 100)
T_homo = T/T_M
#Mu: shear modulus
Mu = Mu_0 * (1 + (T - 300) / T_M * Mu_Temp_Dep)
#sigma_s
sigma_s = np.linspace(math.pow(10, -8), math.pow(10, -3), 100)
sigma_s_homo = sigma_s / Mu
T_homo,sigma_s_homo = np.meshgrid(T_homo, sigma_s_homo)
#D_v
D_v = D_0v * np.exp(-Qv/(R * T))
#D_b
D_b = D_0b * np.exp(-Qb/(R * T))
#D_c
D_c = D_0c * np.exp(-Qc/(R * T))
#D_eff
D_eff = D_v * (10 * Ac / math.pow(b, 2) * np.power(sigma_s_homo / Mu, 2) * D_c / D_v)

#high temperature power-law creep
rate_highTcreep = A2 * D_eff * Mu * b / (R * T)*(np.power(sigma_s_homo / Mu, n))
#low temperature power-law creep
rate_lowTcreep = A2 * D_eff * Mu * b / (R * T)*(np.power(sigma_s_homo, n+2)/(np.power(Mu, n)))
#boundary diffusional flow
rate_highTdiff = 42 * sigma_s_homo * Omega * 3.14 * delta / (R * T * math.pow(d, 3)) * D_b
#lattice diffusional flow
rate_lowTdiff = 42 * sigma_s_homo * Omega / (R * T * math.pow(d, 2))*D_v

'''
#plot rate vs T
plt.figure()
plt.subplot(2,2,1)
plt.plot(T, rate_highTcreep)
plt.xlabel('T')
plt.ylabel('rate_highTcreep')
plt.subplot(2,2,2)
plt.plot(T, rate_lowTcreep)
plt.xlabel('T')
plt.ylabel('rate_lowTcreep')
plt.subplot(2,2,3)
plt.plot(T, rate_highTdiff)
plt.xlabel('T')
plt.ylabel('rate_highTdiff')
plt.subplot(2,2,4)
plt.plot(T, rate_lowTdiff)
plt.xlabel('T')
plt.ylabel('rate_lowTdiff')

#plot rate vs sigma_s
plt.figure()
plt.subplot(2,2,1)
plt.plot(sigma_s, rate_highTcreep)
plt.xlabel('Sigma_s')
plt.ylabel('rate_highTcreep')
plt.subplot(2,2,2)
plt.plot(sigma_s, rate_lowTcreep)
plt.xlabel('Sigma_s')
plt.ylabel('rate_lowTcreep')
plt.subplot(2,2,3)
plt.plot(sigma_s, rate_highTdiff)
plt.xlabel('Sigma_s')
plt.ylabel('rate_highTdiff')
plt.subplot(2,2,4)
plt.plot(sigma_s, rate_lowTdiff)
plt.xlabel('Sigma_s')
plt.ylabel('rate_lowTdiff')
plt.show()
'''

#find the max creep rate
creep_mech=np.empty((4, 100, 100))
mech=np.ones((100, 100))
rate=np.ones((100, 100))
for i in range(100):
    for j in range(100):
        creep_mech[0][i][j] = rate_highTcreep[i][j]
        creep_mech[1][i][j] = rate_lowTcreep[i][j]
        creep_mech[2][i][j] = rate_highTdiff[i][j]
        creep_mech[3][i][j] = rate_lowTdiff[i][j]
        rate[i][j]=np.max([creep_mech[0][i][j],creep_mech[1][i][j],
                              creep_mech[2][i][j],creep_mech[3][i][j]])
        mech[i][j]=np.argmax([creep_mech[0][i][j],creep_mech[1][i][j],
                              creep_mech[2][i][j],creep_mech[3][i][j]])

#select rate
selected_rate = 10**-40
mark = np.empty((100, 100))
for i in range(100):
    for j in range(100):
        if abs(rate[i][j]-selected_rate) < 0.55*selected_rate:
            mark[i][j] = 1
mark = np.ma.masked_where(mark == 0, mark)
#select rate
selected_rate = 10**-25
mark2 = np.empty((100, 100))
for i in range(100):
    for j in range(100):
        if abs(rate[i][j]-selected_rate) < 0.6*selected_rate:
            mark2[i][j] = 1
mark2 = np.ma.masked_where(mark2 == 0, mark2)

#plot
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
ax.set_yscale('log')
#ax.axis('equal')
ax.imshow(mech, cmap='Pastel1',
          extent=[T_homo.min(), T_homo.max(), sigma_s_homo.max(), sigma_s_homo.min()], aspect=1/5)
ax.imshow(mark, cmap='gray',
          extent=[T_homo.min(), T_homo.max(), sigma_s_homo.max(), sigma_s_homo.min()], aspect=1/5)
ax.imshow(mark2, cmap='gray',
          extent=[T_homo.min(), T_homo.max(), sigma_s_homo.max(), sigma_s_homo.min()], aspect=1/5)
ax.invert_yaxis()
ax.set_xlabel('T/Tm')
ax.set_ylabel('Sigma/mu')
fig.tight_layout()
plt.show()


