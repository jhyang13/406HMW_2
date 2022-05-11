import matplotlib.pyplot as plt
import math

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
k = math.pow(10, -23) * 1.38
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
T = [i for i in range(654, 3272)]
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
    k = k + 3 * math.pow(10, -7)
    if j >= 2618:
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
for i in range(0, 2618):
    D_eff.append(D_v[i] * (1 + 10 * Ac / math.pow(b, 2) * math.pow(sigma_s[i] / Mu[i], 2) * D_c[i] / D_v[i]))

#high temperature power-law creep
rate_highTcreep = []
for i in range(0, 2618):
    rate_highTcreep.append(A2 * D_eff[i] * Mu[i] * b / (k * T[i])*(math.pow(sigma_s[i] / Mu[i], n)))
#low temperature power-law creep
rate_lowTcreep = []
for i in range(0, 2618):
    rate_lowTcreep.append(A2 * D_eff[i] * Mu[i] * b / (k * T[i])*(math.pow(sigma_s[i], n+2)/(math.pow(Mu[i],n))))
#boundary diffusional flow
rate_highTdiff = []
for i in range(0, 2618):
    rate_highTdiff.append(42 * sigma_s[i] * Omega * 3.14 * delta / (k * T[i] * math.pow(d, 3)) * D_b[i])
#lattice diffusional flow
rate_lowTdiff = []
for i in range(0, 2618):
    rate_lowTdiff.append(42 * sigma_s[i] * Omega / (k * T[i] * math.pow(d, 2))*D_v[i])

#plot
'''plt.figure()
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
plt.show()'''

#compare and seek the same value between two strain rates
#rate_highTcreep vs rate_lowTcreep
I1 = []
for i in range(0, 2618):
    if rate_highTcreep[i] == rate_lowTcreep[i]:
        I1.append(i)
#rate_highTcreep vs rate_highTdiff
I2 = []
for i in range(0, 2618):
    if rate_highTcreep[i] == rate_highTdiff[i]:
        I2.append(i)
#rate_highTcreep vs rate_lowTdiff
I3 = []
for i in range(0, 2618):
    if rate_highTcreep[i] == rate_lowTdiff[i]:
        I3.append(i)
#rate_lowTcreep vs rate_highTdiff
I4 = []
for i in range(0, 2618):
    if rate_lowTcreep[i] == rate_highTdiff[i]:
        I4.append(i)
#rate_lowTcreep vs rate_lowTdiff
I5 = []
for i in range(0, 2618):
    if rate_lowTcreep[i] == rate_lowTdiff[i]:
        I5.append(i)
#rate_highTdiff vs rate_lowTdiff
I6 = []
for i in range(0, 2618):
    if rate_highTdiff[i] == rate_lowTdiff[i]:
        I6.append(i)
#find the corresponding sigma_s
sigma_s1 = []
T1 = []
for i in I1:
    sigma_s1.append(sigma_s[i])
    T1.append(T[i])
sigma_s2 = []
T2 = []
for i in I2:
    sigma_s2.append(sigma_s[i])
    T2.append(T[i])
sigma_s3 = []
T3 = []
for i in I3:
    sigma_s3.append(sigma_s[i])
    T3.append(T[i])
sigma_s4 = []
T4 = []
for i in I4:
    sigma_s4.append(sigma_s[i])
    T4.append(T[i])
sigma_s5 = []
T5 = []
for i in I5:
    sigma_s5.append(sigma_s[i])
    T5.append(T[i])
sigma_s6 = []
T6 = []
for i in I6:
    sigma_s6.append(sigma_s[i])
    T6.append(T[i])

#plot HMW_2_2
plt.figure()
plt.plot(T1, sigma_s1)
plt.plot(T2, sigma_s2)
plt.plot(T3, sigma_s3)
plt.plot(T4, sigma_s4)
plt.plot(T5, sigma_s5)
plt.plot(T6, sigma_s6)
plt.show()
