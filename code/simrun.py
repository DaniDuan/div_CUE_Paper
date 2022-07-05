from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

N = 5000 # Number of consumers
Tref = 273.15 + 0 # Reference temperature Kelvin
lf = 0.4 # Leakage
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
k = 0.0000862 # Boltzman constant

rho_R = 0
B_R0 = 1.70; B_R0_var = 0.1* B_R0
Ea_R_mean = 0.67; Ea_R_var = 0.04*Ea_R_mean
cov_xy_R = rho_R * B_R0_var**0.5 * Ea_R_var ** 0.5
mean_R = [B_R0, Ea_R_mean]
cov_R = [[B_R0_var, cov_xy_R], [cov_xy_R, Ea_R_var]]  

rho_U = -0.5
B_U0 = 4.47; B_U0_var = 0.1* B_U0
Ea_U_mean = 0.82; Ea_U_var = (0.04*Ea_U_mean)
cov_xy_U = rho_U * B_U0_var**0.5 * Ea_U_var ** 0.5
mean_U = [B_U0, Ea_U_mean]
cov_U = [[B_U0_var, cov_xy_U], [cov_xy_U, Ea_U_var]]  

B_R, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
B_U, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T

# plt.plot(B_U, Ea_U, 'x')
# plt.plot(B_R, Ea_R, 'x')
# plt.show()

plt.hist(Ea_U, 30, color = "orangered", density = True, alpha = 0.5, label = "Uptake rate")
plt.hist(Ea_R, 30, color = "g", density = True, alpha = 0.5, label = "Respiration rate")
plt.axvline(0.82, color="orangered", linestyle='dashed') # Median
plt.axvline(0.67, color='g', linestyle='dashed') # Median
plt.xlabel("Activation Energy (Ea)")
plt.ylabel("Density")
plt.legend()
plt.show()

##############################################################################################
N = 2 # Number of consumers
M = 2 # Number of resources

# Temperature params
T = 273.15 + 20 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration
# covariance of B and Ea, abs(rho)<=1
rho_R = 0.5 
rho_U = 0.5

# Assembly
ass = 1 # Assembly number, i.e. how many times the system can assemble
t_fin = 50 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant

result_array, rich_series, CUE_out, Ea_CUE_out, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U)
print(np.mean(rich_series))

from matplotlib.lines import Line2D

rich = np.array([len(np.where(result_array[i,0:N])[0]) for i in range(len(result_array))]).reshape(ass,t_fin)
rich_mean = np.mean(rich, axis = 0)
rich_ci = 1.96 * np.std(rich,axis = 0)/(ass**0.5)

t_plot = np.linspace(0,len(result_array),len(result_array))

plt.plot(t_plot, result_array[:,N:N+M], 'b-', linewidth=2.5, label = "Resources")
plt.plot(t_plot, result_array[:,0:N], 'g-', linewidth=2.5, label = "Consumers")
plt.ylabel('Consumer & Resource Concentration')
plt.xlabel('Time')
plt.legend([Line2D([0], [0], color='green', lw=2), Line2D([0], [0], color='blue', lw=2)], ['Consumers', 'Resources'])
plt.show()

#######################################################################################################
from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 50 # Number of consumers
M = 100 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration
# covariance of B and Ea, abs(rho)<=1
rho_R = 0
rho_U = -1

# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
t_fin = 2500 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant

T_c = 6 # How many temperatures to cover (how many cycles to run)

########## Running Model ###########
rich = np.empty((0, ass))

for i in range(T_c):
    T = 273.15 + 5 * i # Temperature
    rich_series = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U)[1]
    rich = np.append(rich, [rich_series.flatten()], axis = 0)
    
rich_mean = np.mean(rich, axis = 1)
rich_ci =  1.96 * np.std(rich,axis = 1)/(T_c**0.5)

T_plot = range(0, 5*T_c, 5)

plt.plot(T_plot, rich_mean)
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, color='b', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
plt. show()
