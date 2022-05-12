import matplotlib.pyplot as plt
import numpy as np

#parameters
d = 100*10**-6
omega = 1.8*10**-29
b = 2.86*10**-10
T_m = 3271
n = 4.2
A = 7.5*10**5
R = 8.314
A2 = A*(np.sqrt(3))**(n+1)
delta = 1*10**(-9)
grid = 300

#variables
T = np.linspace(120, T_m, grid)
T_homo = T/T_m
sigma_s_mu = np.logspace(-6, -3, grid)
T,sigma_s_mu = np.meshgrid(T, sigma_s_mu)

#functions
mu_0=6.12*10**4
T_m=3271
mu= mu_0 * (1 + (-0.42) / T_m * (T-300))

D0_v=1.2*10**-5
Qv=413*10**3
Dv=D0_v*np.exp(-Qv/R/T)

D0_b=5.7*10**-14
Qb=280*10**3
Db=D0_b*np.exp(-Qb/R/T) / delta

ac_D0c=1.0*10**-23
Qc=280*10**3
ac_Dc=ac_D0c*np.exp(-Qc/R/T)

#creep rate
#Deff=Dv*(10/(b**2)*ac_Dc/Dv*(sigma_s_mu)**2)
rate_highTcreep=A2*(Dv)*mu*b/R/T*(sigma_s_mu)**n
rate_lowTcreep=A2*10*ac_Dc*mu*b/R/T*(sigma_s_mu)**(n+2)/b**2
rate_highTdiff= (42*omega*sigma_s_mu*mu)/R/T/(d**2)*(Dv)
rate_lowTdiff= (42*omega*sigma_s_mu*mu)/R/T/(d**3)*(Db*3.14159*delta)

creep_mech=np.empty((4,grid,grid))
mech=np.ones((grid,grid))
rate=np.ones((grid,grid))
for i in range(grid):
    for j in range(grid):
        creep_mech[0][i][j]=rate_highTcreep[i][j]
        creep_mech[1][i][j]=rate_lowTcreep[i][j]
        creep_mech[2][i][j]=rate_highTdiff[i][j]
        creep_mech[3][i][j]=rate_lowTdiff[i][j]
        rate[i][j]=np.max([creep_mech[0][i][j],creep_mech[1][i][j],
                              creep_mech[2][i][j],creep_mech[3][i][j]])
        mech[i][j]=np.argmax([creep_mech[0][i][j],creep_mech[1][i][j],
                              creep_mech[2][i][j],creep_mech[3][i][j]])

#plot
fig, ax = plt.subplots(1, 1, figsize=(4,4))
ax.set_yscale('log')
#ax.axis('equal')
ax.imshow(mech,cmap='plasma',
          extent=[T_homo.min(), T_homo.max(), sigma_s_mu.max(), sigma_s_mu.min()], aspect=1/5)
ax.invert_yaxis()
ax.set_xlabel('T/Tm')
ax.set_ylabel('Sigma/mu')
fig.tight_layout()
plt.show()

