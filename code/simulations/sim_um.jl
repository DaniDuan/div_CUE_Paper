include("../sim_frame.jl")

########################### u & m ###############################
N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
# ρ_t=[-0.1384 -0.1384]; # realistic covariance
Tr=273.15+13; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 26

ρ_t = [-0.9999 -0.9999]
Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (2400, 900));
ax1 = Axis(f[1,1], xlabel = "log(Bᵤ)", ylabel = "Eᵤ", xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f[1,2], xlabel = "Temperature (°C)", ylabel = "Uptake (log)", xlabelsize = 50, ylabelsize = 50)
ax3 = Axis(f[2,1], xlabel = "log(Bₘ)", ylabel = "Eₘ", xlabelsize = 50, ylabelsize = 50)
ax4 = Axis(f[2,2], xlabel = "Temperature (°C)", ylabel = "Respiration (log)", xlabelsize = 50, ylabelsize = 50)
T = range(273.15, 273.15+num_temps-1, length = num_temps)
k = 0.0000862 # Boltzman constant
B0 = [log((0.138/(1 - L - 0.22)) * exp((-0.82/k) * ((1/Tr)-(1/273.15)))/(1 + (0.82/(Ed - 0.82)) * exp(Ed/k * (1/308.15 - 1/Tr)))) log(0.138 *exp((-0.67/k) * ((1/Tr)-(1/273.15)))/(1 + (0.67/(Ed - 0.67)) * exp(Ed/k * (1/311.15 - 1/Tr))))]# Using CUE0 = 0.22, mean growth rate = 0.48
B0_var = 0.17*abs.(B0); E_mean = [0.82 0.67]; E_var =  0.14*abs.(E_mean)
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]
for i in 1:N
    allu = rand(MvNormal(meanv[:,1], cov_u), 1)
    allm = rand(MvNormal(meanv[:,2], cov_m), 1)
    B = [exp.(allu[1,:]) exp.(allm[1,:])]
    E = [allu[2,:] allm[2,:]]
    Tpu = 273.15 .+ rand(Normal(35, 5), 1)
    Tpm = Tpu .+ 3
    Tp = [Tpu Tpm]
    temp_p = log.(B .* exp.((-E./k) .* ((1 ./T) .-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k .* (1 ./Tp .- 1 ./T))))
    scatter!(ax1, log.(B)[1], E[1], color = ("#FA8328", 1), markersize = 20)
    lines!(ax2, Temp_rich, temp_p[:,1], color = ("#FA8328",0.75), linewidth = 1)
    scatter!(ax3, log.(B)[2], E[2], color = ("#015845", 1), markersize = 20)
    lines!(ax4, Temp_rich, temp_p[:,2], color = ("#015845",0.75), linewidth = 1)
end
f
save("../result/U_R_var_-1.png", f) 


