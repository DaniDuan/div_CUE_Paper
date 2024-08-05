include("./sim_frame.jl")

N=100
M=50
### Temp params 
ρ_t= [-0.3500 -0.3500]; # realistic covariance
Tr=273.15+10; Ed=3.5
###################################
# Generate MiCRM parameters
tspan = (0.0, 1.5e8)
x0 = vcat(fill(0.1, N), fill(1, M))
Ci = fill(0.1, N)
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31

ρ_t= [-0.9999 -0.9999]; 
Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (2400, 900));
ax1 = Axis(f[1,1], xlabel = "log(Bᵤ)", ylabel = "Eᵤ", xlabelsize = 50, ylabelsize = 50)
ax2 = Axis(f[1,2], xlabel = "Temperature (°C)", ylabel = "Uptake (log)", xlabelsize = 50, ylabelsize = 50)
ax3 = Axis(f[2,1], xlabel = "log(Bₘ)", ylabel = "Eₘ", xlabelsize = 50, ylabelsize = 50)
ax4 = Axis(f[2,2], xlabel = "Temperature (°C)", ylabel = "Respiration (log)", xlabelsize = 50, ylabelsize = 50)
T = range(273.15, 273.15+num_temps-1, length = num_temps)
k = 0.0000862 # Boltzman constant
B0 = [-0.8116 -1.4954]# Using CUE0 = 0.22, mean growth rate = 0.48
E_mean = [0.8146 0.5741]

T = range(273.15, 273.15+num_temps-1, length = num_temps)
Tr = 273.15
B = 10

Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50)
hidexdecorations!(ax, ticklabels = true, ticks = true, grid = true, label = false)
hideydecorations!(ax, ticklabels = true, ticks = true, grid = true, label = false)
### main 
E = 0.2; Ed = 0.1; Tp = 273.15 .+ 12
temp_p = B .* exp.((-E/k) .* ((1 ./T) .-(1/Tr)))./(1 .+ (E/(Ed - E)) .* exp.(Ed/k .* (1 ./Tp .- 1 ./T)))
lines!(ax, Temp_rich, temp_p, color = ("#6B8EDE",0.8), linewidth = 5)
f
save("../result/demo.png", f) 
