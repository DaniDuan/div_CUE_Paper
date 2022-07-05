import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pylab as plt
# from matplotlib.lines import Line2D
import size_temp_funcs as st
import parameters as par
import model_func as mod


######### Main Code ###########

######## Set up parameters ###########

N = 50 # Number of consumers
M = 100 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration
# covariance of B and Ea, abs(rho)<=1
rho_R = 0
rho_U = 0

# Assembly
ass = 1 # Assembly number, i.e. how many times the system can assemble
t_fin = 2500 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant

##### Intergrate system forward #####

def ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U):
    '''
    Main function for the simulation of resource uptake and growth of microbial communities.
    '''

    # Setted Parameters
    k = 0.0000862 # Boltzman constant

    B_R0 = np.log(1.70); B_R0_var = 0.05* B_R0
    Ea_R_mean = 0.67; Ea_R_var = 0.04*Ea_R_mean
    cov_xy_R = rho_R * B_R0_var**0.5 * Ea_R_var ** 0.5
    mean_R = [B_R0, Ea_R_mean]
    cov_R = [[B_R0_var, cov_xy_R], [cov_xy_R, Ea_R_var]]  

    B_U0 = np.log(4.47); B_U0_var = 0.05* B_U0
    Ea_U_mean = 0.82; Ea_U_var = (0.04*Ea_U_mean)
    cov_xy_U = rho_U * B_U0_var**0.5 * Ea_U_var ** 0.5
    mean_U = [B_U0, Ea_U_mean]
    cov_U = [[B_U0_var, cov_xy_U], [cov_xy_U, Ea_U_var]]  

    ### Creating empty array for storing data ###
    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    rich_series = np.empty((0))
    CUE_out = np.empty((0,N))
    Sr = np.empty((0,N))
    Ea_CUE_out = np.empty((0,N))

    for i in range(ass):

        ### Resetting values for every assembly ###

        x0 = np.concatenate((np.full([N], 0.1), np.full([M], 1))) # Starting concentration for resources and consumers

        # Set up Ea (activation energy) and B0 (normalisation constant) based on Tom Smith's observations
        T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
        T_pk_R = T_pk_U + 3
        B_R_log, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
        B_U_log, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T
        B_R = np.exp(B_R_log); B_U = np.exp(B_U_log)

        p = np.repeat(p_value, M)  # Resource input

        # Set up model
        U, R, l = par.params(N, M, T, k, Tref, T_pk_U, T_pk_R, B_U, B_R ,Ma, Ea_U, Ea_R, np.repeat(Ea_D,N), lf) # Uptake
        l_sum = np.sum(l, axis=1)

        # Integration
        t = np.linspace(0,t_fin-1,t_fin) 
        pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

        pops = solve_ivp(mod.metabolic_model, t_span= [0,t_fin], y0=x0, t_eval = t, args = pars, method = 'BDF') # Integrate
        pops.y = np.transpose(np.round(pops.y, 7))

        ### Storing simulation results ###
        result_array = np.append(result_array, pops.y, axis=0)

        ### Analysis ###

        # Richness
        rich_series = np.append(rich_series, [len(np.where(pops.y[t_fin-1,0:N])[0])]) # Get the consumer concentrations from last timestep and find out survivors

        # CUE
        CUE = (U @ (1 - l_sum) - R)/(np.sum(U, axis = 1))
        CUE_out = np.append(CUE_out, [CUE], axis = 0)
        Ea_CUE_out = np.append(Ea_CUE_out, [B_R*(Ea_U - Ea_R)/(B_U*(1 - lf) - B_R)], axis = 0)

        # S*
        Sr = np.append(Sr, [R/(np.sum(U, axis = 1)*(1-lf))], axis = 0)

    return result_array, rich_series, CUE_out, Ea_CUE_out, Sr

result_array, rich_series, CUE_out, Ea_CUE_out, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U)

rich_series