from re import S
from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

N = 5000 # Number of consumers
Tref = 273.15 + 10 # Reference temperature Kelvin
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
import scipy as sc
from scipy import stats

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 10 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration
# covariance of B and Ea, abs(rho)<=1
rho_R = -0.5
rho_U = -0.5

# Assembly
ass = 50 # Assembly times at each temperature
t_fin = 4000 # Number of time steps for each temperature
typ = 1 # Functional response, Type I or II
K = 0.5 # Half saturation constant for Monod equation(Type II)
T_c = 26 # How many temperatures to cover (how many cycles to run)

rich = np.empty((0, ass))
eq = np.empty((0, ass))
eq_sur = np.empty((0, ass))
eq_sur_ci = np.empty((0, ass))
dEa = np.empty((0, ass))
dEa_sur = np.empty((0, ass))
dEa_sur_ci = np.empty((0, ass))
sU = []
sR = []
eU = []
eR = []
U_var = []
extU_var = []
sur_var = np.empty((0))
sur_var_ci = np.empty((0))
ext_var = np.empty((0))
ext_var_ci = np.empty((0))
all_Ea = np.empty((0))
all_Ea_ci = np.empty((0))
sur_CUE = []
sur_Ea = []
sur_eq = []
ext_CUE = []
ext_Ea = []
sur_Sr = []
all_Sr = np.empty((0))
sur_overlap = []
ext_overlap = []
sur_crossf = []
ext_crossf = []

all_U = []
all_R = []
all_CUE = []

for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U)
    rich = np.append(rich, [rich_series.flatten()], axis = 0)
    
    sur = [np.where(result_array[(i+1)*t_fin-1, 0:N]) for i in range(ass)]
    ext = [np.where(result_array[(i+1)*t_fin-1, 0:N] == 0) for i in range(ass)]
    
    sU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    sR.append(np.concatenate([R_out[i][sur[i]] for i in range(len(sur))]).ravel())
    eU.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())
    eR.append(np.concatenate([R_out[i][ext[i]] for i in range(len(ext))]).ravel())

    all_U.append(np.concatenate([np.sum(U_out_total, axis = 1).reshape(ass, N)[i] for i in range(len(sur))]).ravel())
    all_R.append(np.concatenate([R_out[i] for i in range(len(ext))]).ravel())

    U_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][sur[i]] for i in range(len(sur))]).ravel())
    extU_var.append(np.concatenate([np.var(U_out_total,axis = 1).reshape(ass, N)[i][ext[i]] for i in range(len(ext))]).ravel())

    all_CUE.append([CUE_out[i] for i in range(T_c)])
    
    sur_CUE.append(np.concatenate([CUE_out[i][sur[i]] for i in range(len(sur))]))
    sur_var = np.append(sur_var, [np.nanmean([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])])
    sur_var_ci = np.append(sur_var_ci, [1.96 * np.nanstd([np.var(CUE_out[i][sur[i]]) for i in range(len(sur))])/(np.sum(rich_series)**0.5)])
    sur_Ea.append(np.concatenate([Ea_CUE_out[i][sur[i]] for i in range(len(sur))]))
    ext_CUE.append(np.concatenate([CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_Ea.append(np.concatenate([Ea_CUE_out[i][ext[i]] for i in range(len(ext))]))
    ext_var = np.append(ext_var, [np.mean([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])])
    ext_var_ci = np.append(ext_var_ci, [1.96 * np.std([np.var(CUE_out[i][ext[i]]) for i in range(len(ext))])/(np.sum(N-rich_series)**0.5)])
    
    sur_Sr.append(np.concatenate([Sr[i][sur[i]] for i in range(len(sur))]))
    sur_eq.append(np.concatenate([Sr[i][sur[i]] - np.mean(Sr,axis=1)[i] for i in range(len(sur))]))
    all_Sr = np.append(all_Sr, [np.mean(Sr)])
    
    eq_sur = np.append(eq_sur, [np.array([np.mean(Sr[i][sur[i]] - np.mean(Sr,axis=1)[i]) for i in range (ass)])], axis = 0)
    eq_sur_ci = np.append(eq_sur_ci, [np.array([1.96*np.std(Sr[i][sur[i]] - np.mean(Sr,axis=1)[i])/(len(sur[i])**0.5) for i in range (ass)])],axis =0)
    eq = np.append(eq, [np.mean([Sr[i,:] - np.mean(Sr,axis=1)[i] for i in range (ass)], axis = 1)], axis= 0)

    all_Ea = np.append(all_Ea, [np.mean(Ea_CUE_out)])
    all_Ea_ci = np.append(all_Ea_ci, [1.96 * np.std(Ea_CUE_out)/((ass*N)**0.5)])
    dEa_sur = np.append(dEa_sur, [np.array([np.mean(np.abs(Ea_CUE_out[i][sur[i]] - np.mean(Ea_CUE_out,axis=1)[i])) for i in range (ass)])], axis = 0)
    dEa_sur_ci = np.append(dEa_sur_ci, [np.array([1.96*np.std(np.abs(Ea_CUE_out[i][sur[i]] - np.mean(Ea_CUE_out,axis=1)[i]))/(len(sur[i])**0.5) for i in range (ass)])],axis =0)
    dEa = np.append(dEa, [np.mean([np.abs(Ea_CUE_out[i,:] - np.mean(Ea_CUE_out,axis=1)[i]) for i in range (ass)], axis = 1)], axis= 0)
    
    sur_overlap.append(np.concatenate([overlap[i][sur[i]] for i in range(len(sur))]))
    ext_overlap.append(np.concatenate([overlap[i][ext[i]] for i in range(len(ext))]))
    
    sur_crossf.append(np.concatenate([crossf[i][sur[i]] for i in range(len(sur))]))
    ext_crossf.append(np.concatenate([crossf[i][ext[i]] for i in range(len(ext))]))
    

# temp_rich = {'rich':rich, 'eq':eq, 'eq_sur':eq_sur, 'eq_sur_ci': eq_sur_ci, 'dEa': dEa, 'dEa_sur': dEa_sur, 'dEa_sur_ci': dEa_sur_ci, \
#      'sU': sU, 'sR': sR, 'eU': eU, 'eR':eR, 'U_var':U_var, 'extU_var':extU_var, 'sur_var':sur_var, 'sur_var_ci':sur_var_ci,\
#      'ext_var': ext_var, 'ext_var_ci':ext_var_ci, 'all_Ea': all_Ea, 'all_Ea_ci': all_Ea_ci, 'sur_CUE':sur_CUE, \
#      'sur_Ea':sur_Ea, 'ext_CUE':ext_CUE, 'ext_Ea':ext_Ea, 'sur_Sr': sur_Sr, 'all_Sr':all_Sr, 'sur_overlap':sur_overlap, \
#      'ext_overlap':ext_overlap, 'sur_crossf':sur_crossf, 'ext_crossf':ext_crossf}
# np.save('../data/temp_rich.npy', temp_rich) 

# temp_rich = np.load('../data/temp_rich.npy',allow_pickle='TRUE').item()

rich_mean = np.nanmean(rich, axis = 1)
rich_ci =  1.96 * np.nanstd(rich,axis = 1)/(ass**0.5)

CUE_mean = np.array([np.mean(sur_CUE[i]) for i in range(T_c)])
CUE_ci = np.array([1.96 * np.std(sur_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])
CUE_ext_mean = np.array([np.mean(ext_CUE[i]) for i in range(T_c)])
CUE_ext_ci = np.array([1.96 * np.std(ext_CUE[i])/(len(ext_CUE[i])**0.5) for i in range(T_c)])

Ea_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
Ea_ci = np.array([1.96 * np.std(sur_Ea[i])/(len(sur_Ea[i])**0.5) for i in range(T_c)])

Sr_mean = np.array([np.mean(sur_Sr[i]) for i in range(T_c)])
Sr_ci = np.array([1.96 * np.std(sur_Sr[i])/(len(sur_Sr[i])**0.5) for i in range(T_c)])

sU_mean = np.array([np.mean(sU[i]) for i in range(T_c)])
sU_ci = np.array([1.96 * np.std(sU[i])/(len(sU[i])**0.5) for i in range(T_c)])
sR_mean = np.array([np.mean(sR[i]) for i in range(T_c)])
sR_ci = np.array([1.96 * np.std(sR[i])/(len(sR[i])**0.5) for i in range(T_c)])

U_var_mean = np.array([np.mean(U_var[i]) for i in range(T_c)])
U_var_ci = np.array([1.96 * np.std(U_var[i])/(len(U_var[i])**0.5) for i in range(T_c)])
extU_var_mean = np.array([np.mean(extU_var[i]) for i in range(T_c)])
extU_var_ci = np.array([1.96 * np.std(extU_var[i])/(len(extU_var[i])**0.5) for i in range(T_c)])

eU_mean = np.array([np.mean(eU[i]) for i in range(T_c)])
eU_ci = np.array([1.96 * np.std(eU[i])/(len(eU[i])**0.5) for i in range(T_c)])
eR_mean = np.array([np.mean(eR[i]) for i in range(T_c)])
eR_ci = np.array([1.96 * np.std(eR[i])/(len(eR[i])**0.5) for i in range(T_c)])

rich_sur = [[np.repeat(rich[i][j], rich[i][j]) for j in range(ass)] for i in range(T_c)]
rich_ext = [[np.repeat(rich[i][j], N - rich[i][j]) for j in range(ass)] for i in range(T_c)]

overlap_sur_mean = np.array([np.mean(sur_overlap[i]) for i in range(T_c)])
overlap_sur_ci = np.array([1.96 * np.std(sur_overlap[i])/(len(sur_overlap[i])**0.5) for i in range(T_c)])
overlap_ext_mean = np.array([np.mean(ext_overlap[i]) for i in range(T_c)])
overlap_ext_ci = np.array([1.96 * np.std(ext_overlap[i])/(len(ext_overlap[i])**0.5) for i in range(T_c)])

crossf_sur_mean = np.array([np.mean(sur_crossf[i]) for i in range(T_c)])
crossf_sur_ci = np.array([1.96 * np.std(sur_overlap[i])/(len(sur_overlap[i])**0.5) for i in range(T_c)])
crossf_ext_mean = np.array([np.mean(ext_crossf[i]) for i in range(T_c)])
crossf_ext_ci = np.array([1.96 * np.std(ext_overlap[i])/(len(ext_overlap[i])**0.5) for i in range(T_c)])

surEa_mean = np.array([np.mean(sur_Ea[i]) for i in range(T_c)])
surEa_ci = np.array([1.96 * np.std(sur_Ea[i])/(len(sur_Ea[i])**0.5) for i in range(T_c)])

Ea = []
Ea_e = []
for i in range(T_c): 
    n = 0
    for j in rich[i]:
        j = int(j)
        Ea.append(sur_Ea[i][n:n+j])
        Ea_e.append(ext_Ea[i][n:n+j])
        n = n + j

A = list(zip(rich.flatten(),Ea))
B = list(zip(rich.flatten(),Ea_e))
A.sort(key = lambda x: x[0])
B.sort(key = lambda x: x[0])

rich_v = []
for x in rich.flatten():
    if x not in rich_v:
        rich_v.append(x)
rich_v = np.sort(rich_v)


Ea_sorted = []
Ea_e_sorted = []
for i in rich_v:
    sorting = np.empty((0))
    sorting_e = np.empty((0))
    for j in range(len(A)):
        if [x[0] for x in A][j] == i:
            sorting = np.append(sorting, A[j][1])
            sorting_e = np.append(sorting_e, B[j][1])
    Ea_sorted.append(sorting)
    Ea_e_sorted.append(sorting_e)

    
meanEa = []
ciEa = []
meanEa_e = []
ciEa_e = []
for i in range(len(Ea_sorted)):
    meanEa.append(np.mean(Ea_sorted[i]))
    ciEa.append(1.96 * np.std(Ea_sorted[i])/(len(Ea_sorted[i])**0.5))
    meanEa_e.append(np.mean(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)]))
    ciEa_e.append(1.96 * np.std(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)])/(len(Ea_e_sorted[i][np.where(Ea_e_sorted[i]>-500)])**0.5))

T_plot = range(0, T_c, 1)
T_sur = [[np.repeat(T_plot[i], rich[i][j]) for j in range(ass)] for i in range(T_c)]

plt.rcParams["figure.figsize"] = (15,9)
# plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
plt.rcParams.update({'font.size': 25})
plt.rc('xtick', labelsize=25) 
plt.rc('ytick', labelsize=25) 
# plt.rcParams.update(mpl.rcParamsDefault)

plt.plot(T_plot, rich_mean)
plt.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, alpha=0.1, linewidth=2.5)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Richness')
# plt.text(-5,14,'A',fontsize= 'x-large')
# plt.savefig('../result/rich_temp.png')
plt. show()


CUE_var = np.array([np.var(all_CUE[i]) for i in range(T_c)])
# R_var = np.array([np.var(np.log(all_R[i])) for i in range(T_c)])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T_plot, rich_mean,'C0', label = "Richness",linewidth=2.5)
ax1.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, color='C0', alpha=.1)
ln2 = ax2.plot(T_plot, np.log(CUE_var), 'maroon', label = "CUE Variance",linewidth=2.5)
# ln3 = ax2.plot(T_plot, R_var, 'darkseagreen', label = "respiration variance", alpha=.5,linewidth=2.5)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Richness')
ax2.set_ylabel('CUE Variance (log)')
lns = ln1+ln2
ax1.legend(lns, [i.get_label() for i in lns], loc = 2)
plt.tight_layout()
# text(-5,400,'B',fontsize= 'x-large')
# ax1.text(-5,900,'A',fontsize= 'x-large')
plt.savefig('../result/rich_varCUE_temp.png')
plt. show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T_plot, sU_mean,'darkorange', label = "Survivor Uptake Rate",linewidth=2.5)
ax1.fill_between(T_plot, sU_mean - sU_ci, sU_mean + sU_ci, color='darkorange', alpha=.1)
ln3 = ax1.plot(T_plot, eU_mean,'tan', label = "Extinct Uptake Rate", alpha=.7,linewidth=2.5)
ax1.fill_between(T_plot, eU_mean - eU_ci, eU_mean + eU_ci, color='tan', alpha=.3)
ln2 = ax2.plot(T_plot, sR_mean, 'g', label = "Survivor Respiration Rate",linewidth=2.5)
ax2.fill_between(T_plot, sR_mean - sR_ci, sR_mean + sR_ci, color='g', alpha=.1)
ln4 = ax2.plot(T_plot, eR_mean, 'darkseagreen', label = "Extinct Respiration Rate", alpha=.5,linewidth=2.5)
ax2.fill_between(T_plot, eR_mean - eR_ci, eR_mean + eR_ci, color='darkseagreen', alpha=.1)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Uptake Rate')
ax2.set_ylabel('Repiration Rate')
lns = ln1+ln2+ln3+ln4
ax1.legend(lns, [i.get_label() for i in lns], loc = 2)
plt.tight_layout()
# text(-5,400,'B',fontsize= 'x-large')
# ax1.text(-5,900,'A',fontsize= 'x-large')
plt.savefig('../result/U+R_temp.png')
plt. show()

plt.plot(T_plot, CUE_mean, 'r', label = "Survivor",linewidth=2.5)
plt.fill_between(T_plot, CUE_mean - CUE_ci, CUE_mean + CUE_ci, color='r', alpha=.2)
plt.plot(T_plot, CUE_ext_mean, 'dimgrey', label = "Extinct",linewidth=2.5)
plt.fill_between(T_plot,  CUE_ext_mean - CUE_ext_ci,  CUE_ext_mean + CUE_ext_ci, color='dimgrey', alpha=.2)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('CUE')
plt.legend()
# plt.text(-5,0.6,'B',fontsize= 'x-large')
plt.savefig('../result/CUE_temp.png')
plt. show()


plt.scatter(meanEa, rich_v, color = 'b', alpha = 0.7, s = 50, label = 'Survivors')
plt.errorbar(meanEa, rich_v, xerr=ciEa, fmt=',', color = 'b', alpha = 0.7)
plt.scatter(meanEa_e, rich_v, color = 'grey', alpha = 0.7, s = 50, label = 'Extinct')
plt.errorbar(meanEa_e, rich_v, xerr=ciEa_e, fmt=',', color = 'grey', alpha = 0.7)

n = 0
for i in range(T_c):
    plt.scatter(ext_Ea[i], np.concatenate(rich_ext[i]), color = 'lightgrey', alpha = 0.5, s = 50)
    if n == 0: 
        plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s =50, label = 'Survivors')
        n = 1
    else:
        plt.scatter(sur_Ea[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s = 50)
m, b = np.polyfit(np.delete(meanEa,np.where(np.isnan(meanEa))), np.delete(rich_v,np.where(np.isnan(meanEa))), 1)
x = np.arange((np.max(rich_v)-b)/m, (np.min(rich_v)-b)/m, 0.01)
plt.plot(x, m*x + b, color = 'k',linewidth=3)
# plt.text(0.7, 12.5, '$R^2 = $%s' %np.round(r_value**2, 3), fontsize = 22)
plt.xlabel('Thermal Sensitivity of CUE')
plt.ylabel('Richness')
plt.legend()
# plt.text(-0.35,20,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/EaCUE_richness.png')
plt.show()


# plt.plot(T_plot, overlap_sur_mean, color = 'b', alpha = 0.7, label = 'Resource Overlap - Survivor')
# plt.plot(T_plot, overlap_ext_mean, color = 'lightsteelblue', label = 'Resource Overlap - Extinct')
# plt.fill_between(T_plot, overlap_sur_mean - overlap_sur_ci, overlap_sur_mean + overlap_sur_ci, color = 'b', alpha=.1)
# plt.fill_between(T_plot, overlap_ext_mean - overlap_ext_ci, overlap_ext_mean + overlap_ext_ci, color = 'lightsteelblue', alpha=.1)
# plt.plot(T_plot, crossf_sur_mean, color = 'r', alpha = 0.7, label = 'Facilitation - Survivor')
# plt.plot(T_plot, crossf_ext_mean, color = 'rosybrown', label = 'Facilitation - Extinct')
# plt.fill_between(T_plot, crossf_sur_mean - crossf_sur_ci, crossf_sur_mean + crossf_sur_ci, color = 'r', alpha=.1)
# plt.fill_between(T_plot, crossf_ext_mean - crossf_ext_ci, crossf_ext_mean + crossf_ext_ci, color = 'rosybrown', alpha=.1)
# plt.xlabel('Temperature ($^\circ$C)')
# plt.ylabel('Pair-wise Interactions')
# plt.legend(fontsize = 'x-small', framealpha = 0.4)
# # plt.tight_layout()
# # plt.text(-6,0.4,'A',fontsize= 'x-large')
# plt.savefig('../result/comp_coop.png')
# plt.show()

plt.scatter(crossf_sur_mean, overlap_sur_mean, color = "b", alpha = 0.7)
plt.xlabel("facilitation")
plt.ylabel("overlap")
plt.show()

all_cross = np.concatenate(sur_crossf, axis = 0)
all_over = np.concatenate(sur_overlap, axis = 0)
plt.scatter(all_cross, all_over, color = "k", alpha = 0.5)
plt.xlabel("Facilitation")
plt.ylabel("Overlap")
# plt.savefig('../result/scatter_comp_coop')
plt.show()

plt.hist(all_over, color = "b", label = "Resource Overlap", alpha = 0.7)
plt.hist(all_cross, color = "r", label = "Facilitation", alpha = 0.7)
plt.legend()
# plt.savefig('../result/variation_comp_coop')
plt.show()


plt.plot(T_plot, sur_var,'orangered', label = "Survivor")
plt.fill_between(T_plot, sur_var - sur_var_ci, sur_var + sur_var_ci, color='orangered', alpha=.1)
plt.plot(T_plot, ext_var,'grey', label = "Extinct")
plt.fill_between(T_plot, ext_var - ext_var_ci, ext_var + ext_var_ci, color='grey', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Variance of CUE')
plt.legend()
# plt.text(-4.5,0.31,'F',fontsize= 'x-large')
plt.show()
print(sur_var)

plt.plot(T_plot, np.log(U_var_mean),'darkorange', label = "Survivor",linewidth=2.5)
plt.fill_between(T_plot, np.log(U_var_mean-U_var_ci), np.log(U_var_mean+U_var_ci), color='darkorange', alpha=.2)
plt.plot(T_plot, np.log(extU_var_mean),'grey', label = "Extinct")
plt.fill_between(T_plot, np.log(extU_var_mean-extU_var_ci), np.log(extU_var_mean+extU_var_ci), color='grey', alpha=.2)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('Variance of Uptake (log)')
plt.legend(loc = 2,fontsize = 'x-small')
# plt.tight_layout()
# plt.text(-6,7.2,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/VarU.png')
plt.show()

plt.scatter(Sr_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(Sr_mean, rich_mean, xerr=Sr_ci, fmt=',', color = 'k', alpha = 0.7)
plt.xlabel('Community average S*')
plt.ylabel('Richness')
plt.show()


eq_mean = np.nanmean(eq_sur, axis = 1)
eq_ci =  1.96 * np.nanstd(eq_sur,axis = 1)/(ass**0.5)


plt.plot(T_plot, np.abs(eq_mean), 'r', label = 'Competitive exclusion', alpha = 0.7)
plt.fill_between(T_plot, np.abs(eq_mean) - eq_ci, np.abs(eq_mean) + eq_ci, color='r', alpha=.1)
# plt.plot(T_plot, (1 - overlap_sur_mean)/99, 'r', label = 'Resource partitioning', alpha = 0.7)
# plt.fill_between(T_plot, (1 - overlap_sur_mean)/99 - overlap_sur_ci, (1 - overlap_sur_mean)/99 + overlap_sur_ci, color = 'r', alpha=.1)
plt.xlabel('Temperature ($^\circ$C)')
plt.ylabel('|Survivor $S_i^* - \overline{S^*}$|')
# plt.legend()
plt.savefig('../result/S*.png')
plt.show()

# plt.scatter(meanEa, rich_v, color = 'b', alpha = 0.7, s = 50, label = 'Survivors')
# plt.errorbar(meanEa, rich_v, xerr=ciEa, fmt=',', color = 'b', alpha = 0.7)
# plt.scatter(meanEa_e, rich_v, color = 'grey', alpha = 0.7, s = 50, label = 'Extinct')
# plt.errorbar(meanEa_e, rich_v, xerr=ciEa_e, fmt=',', color = 'grey', alpha = 0.7)

n = 0
for i in range(T_c):
    if n == 0: 
        plt.scatter(sur_Sr[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s =50, label = 'Survivors')
        n = 1
    else:
        plt.scatter(sur_Sr[i], np.concatenate(rich_sur[i]), color = 'b', alpha = 0.5, s = 50)
# m, b = np.polyfit(np.delete(meanEa,np.where(np.isnan(meanEa))), np.delete(rich_v,np.where(np.isnan(meanEa))), 1)
# x = np.arange((np.max(rich_v)-b)/m, (np.min(rich_v)-b)/m, 0.01)
# plt.plot(x, m*x + b, color = 'k',linewidth=3)
# plt.text(0.7, 12.5, '$R^2 = $%s' %np.round(r_value**2, 3), fontsize = 22)
plt.xlabel('$S*$')
plt.ylabel('Richness')
plt.legend()
# plt.text(-0.35,20,'B',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/EaCUE_richness.png')
plt.show()


dEa_mean = np.nanmean(dEa_sur, axis = 1)
dEa_ci = 1.96 * np.nanstd(dEa_sur,axis = 1)/(ass**0.5)

m, b, r_value, p_value, std_err = sc.stats.linregress(dEa_mean, rich_mean)
print(r_value**2)
x = np.arange(np.min(dEa_mean), np.max(dEa_mean), 0.001)
plt.plot(x, m*x + b, color = 'k',linewidth=1, alpha = 0.7)
plt.scatter(dEa_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(dEa_mean, rich_mean, xerr=dEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(dEa_mean, rich_mean, yerr=rich_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Ea Difference')
plt.ylabel('Richness')
plt.show()

m, b, r_value, p_value, std_err = sc.stats.linregress(surEa_mean, rich_mean)
print(r_value**2)
x = np.arange(np.min(surEa_mean), np.max(surEa_mean), 0.001)
plt.plot(x, m*x + b, color = 'k',linewidth=1, alpha = 0.7)
plt.scatter(surEa_mean, rich_mean, s=20, color = 'k', alpha = 0.7)
plt.errorbar(surEa_mean, rich_mean, xerr=surEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(surEa_mean, rich_mean, yerr=rich_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Average survivor $E_{a_{CUE}}$')
plt.ylabel('Richness')
plt.text(1, 13, '$R^2 = $%s' %np.round(r_value**2, 3))
plt.show()

# fig, ax1 = plt.subplots()
# ax2 = ax1.twiny()
# ax1.scatter(dEa_sur, rich_mean, s=20, color = 'b', alpha = 0.7, label = 'Ea difference')
# ax1.set_xlabel('Ea Difference')
# ax2.scatter(eq_sur, rich_mean, s=20, color = 'k', alpha = 0.7, label = 'Fitness difference')
# ax2.set_xlabel('Fitness Difference')
# ax1.set_ylabel('Richness')
# plt.show()

plt.scatter(surEa_mean, Sr_mean, s=100, color = 'k', alpha = 0.7)
plt.errorbar(surEa_mean, Sr_mean, xerr=surEa_ci, fmt=',', color = 'k', alpha = 0.3)
plt.errorbar(surEa_mean, Sr_mean, yerr=Sr_ci, fmt=',', color = 'k', alpha = 0.3)
plt.xlabel('Average survivor $E_{a_{CUE}}$')
plt.ylabel('Average survivor $S^*$')
plt.text(-0.1,0.57,'A',fontsize= 'x-large')
# plt.savefig('../thesis/Figures/Eadiff_fitdiff.png')
plt.show()

sU_var = np.array([np.var(np.log(sU[i])) for i in range(T_c)])
eU_var = np.array([np.var(np.log(eU[i])) for i in range(T_c)])

plt.plot(T_plot, sU_var, label = "Uptake Variance_survivors",linewidth=2.5)
plt.plot(T_plot, eU_var, label = "Uptake Variance_extinct",linewidth=2.5)
plt.legend()
plt.show()


#########################################################################################

import numpy as np
import matplotlib.pyplot as plt

N = 20
k = 0.0000862 # Boltzman constant
Tref = 273.15 + 10 # Reference temperature Kelvin, 0 degrees C
T = 273.15 + np.linspace(0,60,61) # Temperatures
Ea_D = 3.5
lf = 0.4

rho_R = -0.75
rho_U = -0.75

B_R0 = np.log(1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref)))) # Using CUE0 = 0.22, mean growth rate = 0.48
B_R0_var = 0.05* B_R0
Ea_R_mean = 0.67; Ea_R_var = 0.04*Ea_R_mean
cov_xy_R = rho_R * B_R0_var**0.5 * Ea_R_var ** 0.5
mean_R = [B_R0, Ea_R_mean]
cov_R = [[B_R0_var, cov_xy_R], [cov_xy_R, Ea_R_var]]  

B_U0 = np.log((1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))))
B_U0_var = 0.05* B_U0
Ea_U_mean = 0.82; Ea_U_var = (0.04*Ea_U_mean)
cov_xy_U = rho_U * B_U0_var**0.5 * Ea_U_var ** 0.5
mean_U = [B_U0, Ea_U_mean]
cov_U = [[B_U0_var, cov_xy_U], [cov_xy_U, Ea_U_var]]  


np.random.seed(0)
T_pk_U = 273.15 + np.random.normal(35, 5, size = N); T_pk_R = T_pk_U + 3
B_R_log, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
B_U_log, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T


for i in range(N):
    U_Sharpe = B_U[i] * np.exp((-Ea_U[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U[i]/(Ea_D - Ea_U[i])) * np.exp(Ea_D/k * (1/T_pk_U[i] - 1/T))) 
    plt.plot(T - 273.15, np.log(U_Sharpe), color = 'darkorange', linewidth=0.7)

plt.xlabel('Temperature ($^\circ$C)', fontsize = 12) 
plt.ylabel('Uptake Rate', fontsize = 12)
plt.show()

for i in range(N):
    R_Sharpe = B_R[i] * np.exp((-Ea_R[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R[i]/(Ea_D - Ea_R[i])) * np.exp(Ea_D/k * (1/T_pk_R[i] - 1/T))) 
    plt.plot(T - 273.15, np.log(R_Sharpe), 'darkgreen', linewidth=0.7)

plt.xlabel('Temperature ($^\circ$C)', fontsize = 12) 
plt.ylabel('Repiration Rate', fontsize = 12)
plt.show()

#######################################################################################################################

from Bacteria_vector_modular import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

########## Setting Parameters ###########
N = 100 # Number of consumers
M = 50 # Number of resources

# Temperature params
Tref = 273.15 + 10 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

rho_R = -1
rho_U = -1

# covariance of B and Ea, abs(rho)<=1
# Assembly
ass = 30 # Assembly number, i.e. how many times the system can assemble
t_fin = 2500 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant

T_c = 26 # How many temperatures to cover (how many cycles to run)

########## Running Model ###########
rich = np.empty((0, ass))
all_CUE = []
all_Sr = []

for i in range(T_c):
    T = 273.15 + i # Temperature
    result_array, rich_series, l, U_out_total, R_out, CUE_out, Ea_CUE_out, overlap, crossf, Sr = ass_temp_run(t_fin, N, M, T, Tref, Ma, ass, Ea_D, lf, p_value, typ, K, rho_R, rho_U)
    all_CUE.append([CUE_out[i] for i in range(T_c)])
    all_Sr.append([Sr[i] for i in range(T_c)])

    rich = np.append(rich, [rich_series.flatten()], axis = 0)

rich_mean = np.mean(rich, axis = 1)
rich_ci =  1.96 * np.std(rich,axis = 1)/(T_c**0.5)

T_plot = range(0, T_c, 1)

CUE_var = np.array([np.var(all_CUE[i]) for i in range(T_c)])
Sr_var = np.array([np.var(all_Sr[i]) for i in range(T_c)])

plt.rcParams["figure.figsize"] = (10,8)
# plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
plt.rcParams.update({'font.size': 25})
plt.rc('xtick', labelsize=25) 
plt.rc('ytick', labelsize=25) 
# plt.rcParams.update(mpl.rcParamsDefault)
# R_var = np.array([np.var(np.log(all_R[i])) for i in range(T_c)])

plt.plot(T_plot, Sr_var)
plt.show()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ln1 = ax1.plot(T_plot, rich_mean,'C0', label = "Richness",linewidth=2.5)
ax1.fill_between(T_plot, rich_mean - rich_ci, rich_mean + rich_ci, color='C0', alpha=.1)
ln2 = ax2.plot(T_plot, np.log(CUE_var), 'maroon', label = "CUE Variance",linewidth=2.5)
# ln3 = ax2.plot(T_plot, R_var, 'darkseagreen', label = "respiration variance", alpha=.5,linewidth=2.5)
ax1.set_xlabel('Temperature ($^\circ$C)')
ax1.set_ylabel('Richness')
ax2.set_ylabel('CUE Variance (log)')
ax2.set_ylim(-6.5,-0.5)
lns = ln1+ln2
ax1.legend(lns, [i.get_label() for i in lns], loc = 2)
# text(-5,400,'B',fontsize= 'x-large')
# ax1.text(-5,900,'A',fontsize= 'x-large')
plt.tight_layout()
plt.savefig('../result/00.png')
plt. show()

##############################################################################################################
import numpy as np
import matplotlib.pyplot as plt

N = 20
k = 0.0000862 # Boltzman constant
Tref = 273.15 + 10 # Reference temperature Kelvin, 0 degrees C
T = 273.15 + np.linspace(0,25,26) # Temperatures
Ea_D = 3.5
lf = 0.4

rho_R = 0
rho_U = 0

B_R0 = np.log(1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref)))) # Using CUE0 = 0.22, mean growth rate = 0.48
B_R0_var = 0.05* B_R0
Ea_R_mean = 0.67; Ea_R_var = 0.04*Ea_R_mean
cov_xy_R = rho_R * B_R0_var**0.5 * Ea_R_var ** 0.5
mean_R = [B_R0, Ea_R_mean]
cov_R = [[B_R0_var, cov_xy_R], [cov_xy_R, Ea_R_var]]  

B_U0 = np.log((1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))))
B_U0_var = 0.05* B_U0
Ea_U_mean = 0.82; Ea_U_var = (0.04*Ea_U_mean)
cov_xy_U = rho_U * B_U0_var**0.5 * Ea_U_var ** 0.5
mean_U = [B_U0, Ea_U_mean]
cov_U = [[B_U0_var, cov_xy_U], [cov_xy_U, Ea_U_var]]  


np.random.seed(0)
T_pk_U = 273.15 + np.random.normal(35, 5, size = N); T_pk_R = T_pk_U + 3
B_R_log, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
B_U_log, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T
B_R = np.exp(B_R_log); B_U = np.exp(B_U_log)

plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams.update({'font.size': 30})
plt.rc('xtick', labelsize=25) 
plt.rc('ytick', labelsize=25) 

var_U = []
for i in range(N):
    U_Sharpe = B_U[i] * np.exp((-Ea_U[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_U[i]/(Ea_D - Ea_U[i])) * np.exp(Ea_D/k * (1/T_pk_U[i] - 1/T))) 
    var = [var, np.variance(U_Sharpe)]
    plt.plot(T - 273.15, np.log(U_Sharpe), color = 'darkorange', linewidth=0.7)

plt.axvline(25, linewidth = 2.5, color = 'darkblue')
plt.xlabel('Temperature ($^\circ$C)') 
plt.ylabel('Uptake Rate (log)')
plt.tight_layout()
plt.savefig('../result/U_r.png')
plt.show()

for i in range(N):
    R_Sharpe = B_R[i] * np.exp((-Ea_R[i]/k) * ((1/T)-(1/Tref)))/(1 + (Ea_R[i]/(Ea_D - Ea_R[i])) * np.exp(Ea_D/k * (1/T_pk_R[i] - 1/T))) 
    plt.plot(T - 273.15, np.log(R_Sharpe), 'darkgreen', linewidth=0.7)

plt.axvline(25, linewidth = 2.5, color = 'darkblue')
plt.xlabel('Temperature ($^\circ$C)') 
plt.ylabel('Repiration Rate (log)')
plt.tight_layout()
plt.savefig('../result/R_r.png')
plt.show()

###############################################################################################################

N = 100 * 26
k = 0.0000862 # Boltzman constant
Tref = 273.15 + 10 # Reference temperature Kelvin, 0 degrees C
Ea_D = 3.5
lf = 0.4

rho_R = -0.5
rho_U = -0.5

B_R0 = np.log(1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref)))) # Using CUE0 = 0.22, mean growth rate = 0.48
B_R0_var = 0.05* B_R0
Ea_R_mean = 0.67; Ea_R_var = 0.04*Ea_R_mean
cov_xy_R = rho_R * B_R0_var**0.5 * Ea_R_var ** 0.5
mean_R = [B_R0, Ea_R_mean]
cov_R = [[B_R0_var, cov_xy_R], [cov_xy_R, Ea_R_var]]  

B_U0 = np.log((1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref))))
B_U0_var = 0.05* B_U0
Ea_U_mean = 0.82; Ea_U_var = (0.04*Ea_U_mean)
cov_xy_U = rho_U * B_U0_var**0.5 * Ea_U_var ** 0.5
mean_U = [B_U0, Ea_U_mean]
cov_U = [[B_U0_var, cov_xy_U], [cov_xy_U, Ea_U_var]]  


np.random.seed(0)
T_pk_U = 273.15 + np.random.normal(35, 5, size = N); T_pk_R = T_pk_U + 3
B_R_log, Ea_R = np.random.multivariate_normal(mean_R, cov_R, N).T
B_U_log, Ea_U = np.random.multivariate_normal(mean_U, cov_U, N).T
B_R = np.exp(B_R_log); B_U = np.exp(B_U_log)

plt.scatter(B_U, Ea_U)
plt.xlabel("$B_0$")
plt.ylabel("E")
plt.savefig('../result/BU_Ea.png')
plt.show()


########################################################################################
A = 0.52; A_var = 0.001
B = 0.65; B_var = 0.001
cov_xy_AB = 0.6 * A_var**0.5 * B_var ** 0.5
mean_AB = [A, B]
cov_AB = [[A_var, cov_xy_AB], [cov_xy_AB, B_var]]  
all_A, all_B = np.random.multivariate_normal(mean_AB, cov_AB, N).T

# plt.hist(all_B)
# plt.show()

plt.scatter(all_A, all_B)
plt.xlabel("A")
plt.ylabel("B")
plt.savefig("../result/test.png")
plt.show()
