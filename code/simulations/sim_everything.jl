include("../sim_frame.jl")

N=100
M=50
L = 0.3
### Temp params 
# T=15+273.15; 
ρ_t=[-0.1384 -0.1384]; # realistic covariance
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

# ρ_t = [0.0 0.0]
everything = zeros(Float64, num_temps, 16)
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i 
    all = Float64[]; ϵ_sur = Float64[]; ϵ_ext = Float64[]; # Eϵ_sur = Float64[]; Eϵ_ext = Float64[]; 
    ϵ_var = Float64[]; u_sur = Float64[]; u_ext = Float64[]; m_sur = Float64[]; m_ext = Float64[]
    for j in 1:50
        ## generate params
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## Calc CUE
        ϵ = (p.u * x0[N+1:N+M] .* (1-L) .- p.m) ./ (p.u * x0[N+1:N+M])
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        bm = sol.u[length(sol.t)][1:N]
        push!(all, length(bm[bm.>1e-7]))
        append!(ϵ_sur, ϵ[bm.>1e-7]); append!(ϵ_ext, ϵ[bm.<=1e-7]) #CUE
        #append!(Eϵ_sur, p.Eϵ[bm.>1e-7]); append!(Eϵ_ext, p.Eϵ[bm.<=1e-7]) #Ea_CUE
        append!(u_sur, sum(p.u, dims = 2)[bm.>1e-7]); append!(u_ext, sum(p.u, dims = 2)[bm.<=1e-7])
        append!(m_sur, p.m[bm.>1e-7]); append!(m_ext, p.m[bm.<=1e-7])
        push!(ϵ_var, log(var(ϵ)))
    end
    rich = mean(all);rich_err = std(all)/sqrt(length(all))
    everything[Int(i+1),:] = [ mean(all), std(all)/sqrt(length(all)), 
        mean(ϵ_sur), std(ϵ_sur)/sqrt(length(ϵ_sur)), mean(ϵ_ext), std(ϵ_ext)/sqrt(length(ϵ_ext)), mean(ϵ_var), std(ϵ_var)/sqrt(length(ϵ_var)), 
        #mean(Eϵ_sur), std(Eϵ_sur)/sqrt(length(Eϵ_sur)), mean(Eϵ_ext), std(Eϵ_ext)/sqrt(length(Eϵ_ext)),
        mean(u_sur), std(u_sur)/sqrt(length(u_sur)), mean(u_ext), std(u_ext)/sqrt(length(u_ext)), 
        mean(m_sur), std(m_sur)/sqrt(length(m_sur)), mean(m_ext), std(m_ext)/sqrt(length(m_ext))]
    print(i, " °C Complete, ", "richness ",rich,"\n")
end 

col_names = ["richness", "richness_err", "ϵ_sur_mean", "ϵ_sur_err", "ϵ_ext_mean", "ϵ_ext_err", "ϵ_var_mean", "ϵ_var_err", #"Eϵ_sur_mean", "Eϵ_sur_err", "Eϵ_ext_mean", "Eϵ_ext_err",
            "u_sur_mean", "u_sur_err", "u_ext_mean", "u_ext_err", "m_sur_mean", "m_sur_err", "m_ext_mean", "m_ext_err"]
everything = DataFrame(everything, col_names)
# CSV.write("../data/temp_gradient.csv", everything, writeheader=false)
# everything = CSV.read("../data/temp_gradient.csv", DataFrame, header=false)
# rename!(everything, col_names)

# richness plots
Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness", limits = (nothing, nothing, 2.5, 15.5), xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "CUE Variance (log)", limits = (nothing, nothing, -7.5, 1.2), xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, everything.richness, color = ("#6B8EDE",0.8), linewidth = 5, label = "Richness")
band!(ax1, Temp_rich, everything.richness .- everything.richness_err, everything.richness .+ everything.richness_err, color = ("#6B8EDE", 0.2))
lines!(ax2, Temp_rich, everything.ϵ_var_mean, color = ("#EF8F8C", 0.8), linewidth = 5, label = "CUE Variance")
band!(ax2, Temp_rich,  everything.ϵ_var_mean .- everything.ϵ_var_err, everything.ϵ_var_mean .+ everything.ϵ_var_err, color = ("#EF8F8C", 0.2))
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#6B8EDE",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["Richness", "CUE Variance"], halign = :left, valign = :top)
f

# save("../result/00.png", f) 

# save("../result/rich_temp.png", f) 

f = Figure(fontsize = 35, resolution = (1200, 900));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "CUE", xlabelsize = 50, ylabelsize = 50)
lines!(ax, Temp_rich, everything.ϵ_sur_mean, color = ("#EF8F8C",1), linewidth = 5, label = "Survivor")
band!(ax, Temp_rich, everything.ϵ_sur_mean .- everything.ϵ_sur_err , everything.ϵ_sur_mean .+ everything.ϵ_sur_err , color = ("#EF8F8C", 0.2))
lines!(ax, Temp_rich, everything.ϵ_ext_mean, color = ("#4F363E", 0.6), linewidth = 5, label = "Extinct")
band!(ax, Temp_rich,  everything.ϵ_ext_mean .- everything.ϵ_ext_err, everything.ϵ_ext_mean .+ everything.ϵ_ext_err, color = ("#4F363E", 0.2))
axislegend(position = :rb)
f
# save("../result/CUE_temp.png", f) 

Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Uptake Rate", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false) 
ax2 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Maintenance Respiration", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false) 
lines!(ax1, Temp_rich, everything.u_sur_mean, color = ("#FA8328", 1), linewidth = 5)
band!(ax1, Temp_rich, everything.u_sur_mean .- everything.u_sur_err, everything.u_sur_mean .+ everything.u_sur_err, color = ("#FA8328", 0.2))
lines!(ax1, Temp_rich, everything.u_ext_mean, color = ("#FA8328", 0.4), linewidth = 5)
band!(ax1, Temp_rich, everything.u_ext_mean .- everything.u_ext_err, everything.u_ext_mean .+ everything.u_ext_err, color = ("#FA8328", 0.2))
lines!(ax2, Temp_rich, everything.m_sur_mean, color = ("#015845", 1), linewidth = 5)
band!(ax2, Temp_rich,  everything.m_sur_mean .- everything.m_sur_err, everything.m_sur_mean .+ everything.m_sur_err, color = ("#015845", 0.2))
lines!(ax2, Temp_rich, everything.m_ext_mean, color = ("#015845", 0.4), linewidth = 5)
band!(ax2, Temp_rich,  everything.m_ext_mean .- everything.m_ext_err, everything.m_ext_mean .+ everything.m_ext_err, color = ("#015845", 0.2))
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#FA8328", 1), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color =("#015845", 1), linestyle = nothing, linewidth = 5)]
l3 = [LineElement(color = ("#FA8328", 0.4), linestyle = nothing, linewidth = 5)]
l4 = [LineElement(color = ("#015845", 0.4), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2, l3, l4], tellheight = false, tellwidth = false, ["Survivor uptake", "Survivor maintenance", "Extinct uptake", "Extinct maintenance"], halign = :left, valign = :top)
f
# save("../result/U+R_temp.png", f) 
