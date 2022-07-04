from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

N = 1000 # Number of consumers
Tref = 273.15 + 0 # Reference temperature Kelvin
lf = 0.4 # Leakage
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
k = 0.0000862 # Boltzman constant

B_R0 = 1.70 # 1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref))) # Using CUE0 = 0.22, mean growth rate = 0.48
Ea_R_int = 0.67
rho_R = 1 # abs(rho)<=1
cov_xy_R = rho_R * (0.1* B_R0)**0.5 * (0.1*Ea_R_int) ** 0.5
mean_R = [B_R0, Ea_R_int]
cov_R = [[(0.1* B_R0), cov_xy_R], [cov_xy_R, (0.1*Ea_R_int)]]  

B_U0 = 4.47 # (1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref)))
Ea_U_int = 0.82
rho_U = 1 # abs(rho)<=1
cov_xy_U = rho_U * (0.1* B_U0)**0.5 * (0.1*Ea_U_int) ** 0.5
mean_U = [B_U0, Ea_U_int]
cov_U = [[(0.1* B_U0), cov_xy_U], [cov_xy_U, (0.1*Ea_U_int)]]  

B_R, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
B_U, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T

plt.plot(B_U, Ea_U, 'x')
plt.plot(B_R, Ea_R, 'x')
plt.show()