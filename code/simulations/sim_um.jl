include("./sim_frame.jl")

########################### u & m ###############################
N=100
M=50
L = 0.3
### Temp params 
Tr=273.15+10; Ed=3.5
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
x0 = vcat(fill(0.1, N), fill(1, M))
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
num_temps = 31

# ρ_t= [-0.3500 -0.3500]; # realistic covariance
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
B0_var = 0.17 .* abs.(B0); E_mean = [0.8146 0.5741]; E_var =  0.1364 .* E_mean
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]

allNu = zeros(Float64, N, 31)
allNm = zeros(Float64, N, 31)
allB = zeros(Float64, N, 2)
allE = zeros(Float64, N, 2)
for i in 1:N
    allu = rand(MvNormal(meanv[:,1], cov_u), 1)
    allm = rand(MvNormal(meanv[:,2], cov_m), 1)
    B = [exp.(allu[1,:]) exp.(allm[1,:])]
    E = [allu[2,:] allm[2,:]]
    Tpu = 273.15 .+ rand(Normal(35, 5), 1)
    Tpm = Tpu .+ 2.6
    Tp = [Tpu Tpm]
    temp_p = log.(B .* exp.((-E./k) .* ((1 ./T) .-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k .* (1 ./Tp .- 1 ./T))))
    scatter!(ax1, log.(B)[1], E[1], color = ("#FA8328", 1), markersize = 20)
    lines!(ax2, Temp_rich, temp_p[:,1], color = ("#FA8328",0.75), linewidth = 1)
    scatter!(ax3, log.(B)[2], E[2], color = ("#015845", 1), markersize = 20)
    lines!(ax4, Temp_rich, temp_p[:,2], color = ("#015845",0.75), linewidth = 1)
    allNu[i, :] = temp_p[:,1]
    allNm[i, :] = temp_p[:,2]
    allB[i, :] = B
    allE[i, :] = E
end
# argmin(var(allNu, dims = 1))
# argmin(var(allNm, dims = 1))
# vlines!(ax2, Temp_rich[argmin(var(allNu, dims = 1))[2]], color = :red, linestyle = :dash)
# vlines!(ax4, Temp_rich[argmin(var(allNm, dims = 1))[2]], color = :red, linestyle = :dash)


cor(log.(allB[:,1]),allE[:,1])
cor(log.(allB[:,2]),allE[:,2])

f
# save("../result/BEum.png", f) 
save("../result/U_R_var_-1.png", f) 




############### Plotting u+m and CUE #############
ρ_t=[-0.3500 -0.3500]
num_temps = 31
Temp_rich = range(0, num_temps-1, length = num_temps)
T = Temp_rich .+ 273.15
k = 0.0000862 # Boltzman constant
B0 = [log(0.2875 * exp((-0.82/k) * ((1/Tr)-(1/273.15)))/(1 + (0.82/(Ed - 0.82)) * exp(Ed/k * (1/308.15 - 1/Tr)))) log(0.138 *exp((-0.67/k) * ((1/Tr)-(1/273.15)))/(1 + (0.67/(Ed - 0.67)) * exp(Ed/k * (1/311.15 - 1/Tr))))]# Using CUE0 = 0.22, mean growth rate = 0.48
B0_var = 0.17*abs.(B0); E_mean = [0.82 0.67]; E_var =  0.14*abs.(E_mean)
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] co v_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]
allu = rand(MvNormal(meanv[:,1], cov_u), 1)
allm = rand(MvNormal(meanv[:,2], cov_m), 1)
B = [exp.(allu[1,:]) exp.(allm[1,:])]
E = [allu[2,:] allm[2,:]]
Tpu = 273.15 .+ rand(Normal(35, 5), 1)
Tpm = Tpu .+ 3
Tp = [Tpu Tpm]
temp_p = B .* exp.((-E./k) .* ((1 ./T) .-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k .* (1 ./Tp .- 1 ./T)))
ϵ = (temp_p[:,1] .* (1 - L) .- temp_p[:,2]) ./ (temp_p[:,1])

f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Metabolic Rate", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false) 
ax2 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Carbon Use Efficiency", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false) 
lines!(ax1, Temp_rich, temp_p[:,1], color = ("#FA8328",1), linewidth = 5)
lines!(ax1, Temp_rich, temp_p[:,2], color = ("#015845",1), linewidth = 5)
lines!(ax2, Temp_rich, ϵ, color = ("#EF8F8C",1), linewidth = 5)
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#FA8328",1), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color =("#015845",1), linestyle = nothing, linewidth = 5)]
l3 = [LineElement(color =("#EF8F8C",1), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2, l3], tellheight = false, tellwidth = false, ["Resource uptake", "Maintenance respiration", "Carbon use efficiency"], halign = :left, valign = :top)
f
Label(f[1,1, TopLeft()], "(a)")

save("../result/umCUE_in.png", f) 


######### Fitting ϵ to temperature performance curves ############
num_temps = 51
k = 0.0000862 # Boltzman constant
T = range(273.15, 273.15+num_temps-1, length = num_temps)
B0 = [log((0.138/(1 - L - 0.22)) * exp((-0.82/k) * ((1/Tr)-(1/273.15)))/(1 + (0.82/(Ed - 0.82)) * exp(Ed/k * (1/308.15 - 1/Tr)))) log(0.138 *exp((-0.67/k) * ((1/Tr)-(1/273.15)))/(1 + (0.67/(Ed - 0.67)) * exp(Ed/k * (1/311.15 - 1/Tr))))]# Using CUE0 = 0.22, mean growth rate = 0.48
B0_var = 0.17*abs.(B0); E_mean = [0.82 0.67]; E_var =  0.14*abs.(E_mean)
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]

# The SS TPC model
temp_SS(T, params) = params[1] .* exp.((-params[2]./k) * ((1/T)-(1/Tr)))./(1 .+ (params[2]./(params[4] .- params[2])) .* exp.(params[4]/k * (1 ./params[3] .- 1/T)))

ϵ = Matrix{Float64}(undef, N, num_temps)
for i in 1: num_temps
    Random.seed!(0)
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, T=T[i], ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ϵ_T = (p.u * x0[N+1:N+M] .* (1 .- p.L) .- p.m) ./ (p.u * x0[N+1:N+M])
    ϵ[:,i] = ϵ_T
end 

Tpϵ = []
for i in 1:N
    push!(Tpϵ, T[argmax(ϵ[i,:])])
end 


N = 200
ρ_t= [-0.3500 -0.3500] #[-0.1384 -0.1384]
num_temps = 31
Temp_rich = range(0, num_temps-1, length = num_temps)
T = Temp_rich .+ 273.15
k = 0.0000862 # Boltzman constant
B0 = [-0.8116 -1.4954]
B0_var = 0.17 .* abs.(B0); E_mean = [0.8146 0.5741]; E_var =  0.1364 .* E_mean
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]
allu = rand(MvNormal(meanv[:,1], cov_u), N)
allm = rand(MvNormal(meanv[:,2], cov_m), N)
B = [exp.(allu[1,:]) exp.(allm[1,:])]
E = [allu[2,:] allm[2,:]]
cor_value_u = cor(allu[1,:], E[:,1]); println(cor_value_u)
cor_value_m = cor(allm[1,:], E[:,2]); println(cor_value_m)

# 
f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "log(B0)", ylabel = "E", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, allu[1,:], allu[2,:], color = ("#FA8328", 1), markersize = 25, label = "Resource uptake rate")
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(a)")
f
save("../result/Bu_Ea_em.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "log(B0)", ylabel = "E", xlabelsize = 50, ylabelsize = 50)
scatter!(ax, allm[1,:], allm[2,:], color = ("#015845", 1), markersize = 25, label = "Respiration rate")
axislegend(position = :rt)
Label(f[1,1, TopLeft()], "(b)")
f
save("../result/Bm_Ea_em.png", f) 


#### only u and m

ρ_t=[-0.3500 -0.3500]
num_temps = 61
Temp_rich = range(0, num_temps-1, length = num_temps)
T = Temp_rich .+ 273.15
k = 0.0000862 # Boltzman constant
B0 = [log(0.2875 * exp((-0.82/k) * ((1/Tr)-(1/273.15)))/(1 + (0.82/(Ed - 0.82)) * exp(Ed/k * (1/308.15 - 1/Tr)))) log(0.138 *exp((-0.67/k) * ((1/Tr)-(1/273.15)))/(1 + (0.67/(Ed - 0.67)) * exp(Ed/k * (1/311.15 - 1/Tr))))]# Using CUE0 = 0.22, mean growth rate = 0.48
B0_var = 0.17*abs.(B0); E_mean = [0.82 0.67]; E_var =  0.14*abs.(E_mean)
cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
meanv = [B0 ; E_mean]
cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]


f = Figure(fontsize = 35, resolution = (1200, 800));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Uptake Rate", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false) 
for i in 1:N
    allu = rand(MvNormal(meanv[:,1], cov_u), 1)
    allm = rand(MvNormal(meanv[:,2], cov_m), 1)
    B = [exp.(allu[1,:]) exp.(allm[1,:])]
    E = [allu[2,:] allm[2,:]]
    Tpu = 273.15 .+ rand(Normal(35, 5), 1)
    Tpm = Tpu .+ 2.6
    Tp = [Tpu Tpm]
    temp_p = B .* exp.((-E./k) .* ((1 ./T) .-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k .* (1 ./Tp .- 1 ./T)))
    lines!(ax1, Temp_rich, temp_p[:,1], color = ("#FA8328",0.75), linewidth = 1)
    # lines!(ax1, Temp_rich, temp_p[:,2], color = ("#015845",0.75), linewidth = 1)
end
vlines!(ax1, 30, color = :black, linestyle = :dash)
f

save("../result/u.png", f) 

N = 1
M = 5
### Heatmaps 
T = 273.15 + 10
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)

f = Figure(resolution = (1200, 1200));
ax = Axis(f[1,1], ygridvisible = true, xgridvisible = true)

hidedecorations!(ax, ticks = false)
heatmap!(ax, p.l[1, :,:], colormap = :grays)
f

save("../result/li.png", f) 
