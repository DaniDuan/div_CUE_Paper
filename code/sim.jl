# cd("WORKING_DIRECTORY")

# Load libraries
using Distributions
using LinearAlgebra
using DifferentialEquations
# using Plots, StatsPlots
using Sundials
using Parameters
using CSV, DataFrames
using CairoMakie
using LsqFit 


# Include simulation code files
include("micrm_params.jl") # Contains function gereate_params with default sampling scheme

include("dx_v2.jl") # Defines differential equations for MiCRM integration use dxx instead of dx

include("LV_dx.jl") # Defines LV differential equatons, use LV_dx

include("temp.jl")

####################################################################################################################################################################################
function F_m(N, M, kw)
    if haskey(kw, :T)
        m = kw[:tt][:,2]
    else 
        m = fill(0.2, N)
    end 
    return m
end

function F_ρ(N, M, kw)
    ρ = fill(1, M)
    return ρ
end

function F_ω(N, M, kw)
    ω = fill(0.0, M)
    return ω
end

function F_u(N, M, kw)
    if haskey(kw, :T)
        u_sum = kw[:tt][:,1]
    else 
        u_sum = fill(2.5, M) 
    end
    diri = transpose(rand(Dirichlet(ones(M)),N))
    u = diri.*u_sum
    return u
end

################################################################################################################################################################################

N=100
M=50
### Temp params 
# T=15+273.15; 
ρ_t=[-0.1384 -0.1384]; Tr=273.15+13; Ed=3.5 #[-0.1384 -0.1384]
###################################
# Generate MiCRM parameters
tspan = (0.0, 15000.0)
# here we define a callback that terminates integration as soon as system reaches steady state
condition(du, t, integrator) = norm(integrator(t, Val{1})) <= eps()
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

all = Float64[]
richness = Float64[]; richness_err = Float64[]
for i in range(0, stop = 25, length = 26)
    T = 273.15 + i 
    # all = Float64[]
    for j in 1:50
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=0.3, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        x0 = vcat(fill(0.1, N), fill(1, M))
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        bm = sol.u[length(sol.t)][1:N]
        push!(all, length(bm[bm.>1e-7]))
        print(i, " °C Completed ",j*2, "% \n")
    end
    # rich = mean(all);rich_err = std(all)/sqrt(length(all))
    # push!(richness, rich); push!(richness_err, rich_err)
    # print(i, " °C Complete, ", "richness ",rich,"\n")
end 

# richness
Temp = Any[]
[append!(Temp, fill(i, 50)) for i in 0:25][1]
Temp = [float(x) for x in Temp if x isa Number]
# plot_scatter =scatter(Temp, all, color = :forestgreen, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.7)
richness_all = all
relative_all = all./maximum(all)
# all_save = DataFrame(richness = all)
# CSV.write("../data/all.csv", all_save, writeheader=false)
# richness_all = CSV.read("../data/all.csv", DataFrame, header=false)[:,1]
# relative_all = richness_all./maximum(richness_all)
################## reading EMP data ###############################
EMP_data = CSV.read("../data/EMP_filtered.csv", DataFrame, header= true)
EMP_relative_rich = EMP_data[:, "relative_rich"]
EMP_temp = EMP_data[:, "Temp"]

EMP = [EMP_temp EMP_relative_rich]
# scatter!(EMP_temp, EMP_relative_rich, color = :darkorange, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.7, label = "EMP")

# quan = DataFrame()
# for i in unique(EMP_data[:, "Temp"])
#     sub = filter(:Temp => n -> n == i, EMP_data)
#     quantile_99 = quantile(sub.relative_rich, 0.99)
#     filtered_df = filter(row -> row.relative_rich > quantile_99, sub)
#     quan = vcat(quan, filtered_df)
# end 

################### binning by integers of temp ##############
filtered_data = copy(EMP_data)
filtered_data.Temp = round.(Int, filtered_data.Temp)
# bins created on unique integer temps
bins = unique(filtered_data.Temp)
filtered_df = DataFrame()
for unique_bin in bins
    # Filter rows for bin
    bin_data = filter(row -> row.Temp == unique_bin, filtered_data)
    # Calculate 99% quantile of richness for bin
    quantile_97_richness = quantile(bin_data.Richness, 0.97)
    # Filter rows with richness above 99% quantile
    filtered_rows = filter(row -> row.Richness > quantile_97_richness, bin_data)
    append!(filtered_df, filtered_rows)
end
println(filtered_df)
################################## fitting EMP to gaussian #################################
gaussian(x, params) = params[1] .* exp.(-(x .- params[2]).^2 ./ (2 * params[3]^2))

initial_params = [1.0, 10, 1.0]
residual(params) = y_data - gaussian.(x_data, params)

# x_data = unique(EMP_data[:, "Temp"]); y_data = quan
x_data = filtered_df[:, "Temp"]; y_data = filtered_df[:, "relative_rich"]
# fitting gaussian
fit_result = curve_fit(gaussian, x_data, y_data, initial_params)

best_fit_params = coef(fit_result)
std_errors = stderror(fit_result)
# Generate x y values for plot
x_fit = range(minimum(x_data), stop=maximum(x_data), length=100)
y_fit = gaussian(x_fit, best_fit_params)

################# fitting the prediction ###############
filtered_pre = DataFrame(Temp = Temp, Richness = richness_all, relative_richness = relative_all)
# bins created on unique integer temps
bins = unique(filtered_pre.Temp)

filtered_df_pre = DataFrame()
for unique_bin in bins
    # Filter rows for bin
    bin_data = filter(row -> row.Temp == unique_bin, filtered_pre)
    # Calculate 99% quantile of rich for bin
    quantile_97_richness = quantile(bin_data.Richness, 0.97)
    # Filter rows with rich above 99% quantile
    filtered_rows = filter(row -> row.Richness > quantile_97_richness, bin_data)
    append!(filtered_df_pre, filtered_rows)
end

x_data_pre = filtered_df_pre[:, "Temp"]; y_data_pre = filtered_df_pre[:, "relative_richness"]
x_data_pre = [float(x) for x in x_data_pre if x isa Number]

# fitting gaussian
fit_result = curve_fit(gaussian, x_data_pre, y_data_pre, initial_params)

best_fit_params = coef(fit_result)
std_errors = stderror(fit_result)
# Generate x y values for plot
x_fit_pre = range(minimum(x_data_pre), stop=maximum(x_data_pre), length=100)
y_fit_pre = gaussian(x_fit_pre, best_fit_params)


########################################## all plots #####################
f = Figure(fontsize = 30, resolution = (1200, 800));
ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Relative Richness")
scatter!(ax, Temp, relative_all, color = "#9497C5", markersize = 7, label = "Prediction")
scatter!(ax, EMP_temp, EMP_relative_rich, color = "#F4D662", markersize = 7, label = "EMP")
scatter!(ax, x_data, y_data, color = "#F8BA17", markersize = 12, fillalpha = 0.7)
scatter!(ax, x_data_pre, y_data_pre, color = "#3878C1", markersize = 12, fillalpha = 0.7)
# Gaussian curves
lines!(ax, x_fit, y_fit, color = ("#F8BA17", 0.8), linewidth = 5, label = "EMP_Gaussian")
lines!(ax, x_fit_pre, y_fit_pre, color = ("#3878C1",0.8), linewidth = 5, label = "MiCRM_Gaussian")
# axislegend(labelsize=10)
f[1, 2] = Legend(f, ax, framevisible = false, labelsize=20)
f

CairoMakie.activate!(type = "png")
save("../result/MiCRM_vs_EMP.png", f) 


####################################################################
N=100
M=50
l_α = 0.3
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

ρ_t = [0.9999 0.9999]
all = Float64[]; ϵ_sur = Float64[]; ϵ_ext = Float64[]; ϵ_var = Float64[];
everything = zeros(Float64, num_temps, 8)
for i in range(0, stop = num_temps-1, length = num_temps)
    T = 273.15 + i 
    all = Float64[]; ϵ_sur = Float64[]; ϵ_ext = Float64[]; ϵ_var = Float64[]
    for j in 1:50
        ## generate params
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=l_α, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        ## Calc CUE
        ϵ = (p.u * x0[N+1:N+M] .* (1-l_α) .- p.m) ./ (p.u * x0[N+1:N+M])
        ## run simulation
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        bm = sol.u[length(sol.t)][1:N]
        push!(all, length(bm[bm.>1e-7]))
        append!(ϵ_sur, ϵ[bm.>1e-7]); append!(ϵ_ext, ϵ[bm.<=1e-7])
        push!(ϵ_var, log(var(ϵ)))
    end
    rich = mean(all);rich_err = std(all)/sqrt(length(all))
    everything[Int(i+1),:] = [rich, rich_err, mean(ϵ_sur), std(ϵ_sur)/sqrt(length(ϵ_sur)), mean(ϵ_ext), std(ϵ_ext)/sqrt(length(ϵ_ext)), mean(ϵ_var), std(ϵ_var)/sqrt(length(ϵ_var))]
    print(i, " °C Complete, ", "richness ",rich,"\n")
end 

col_names = ["richness", "richness_err", "ϵ_sur_mean", "ϵ_sur_err", "ϵ_ext_mean", "ϵ_ext_err", "ϵ_var_mean", "ϵ_var_err"]
everything = DataFrame(everything, col_names)
# CSV.write("../data/temp_gradient.csv", everything, writeheader=false)
# everything = CSV.read("../data/temp_gradient.csv", DataFrame, header=false)
# rename!(everything, col_names)

# richness plots
Temp_rich = range(0, num_temps-1, length = num_temps)
f = Figure(fontsize = 35, resolution = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness")
ax1 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "Richness", xlabelsize = 50, ylabelsize = 50, ygridvisible = false, xgridvisible = false)
ax2 = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "CUE Variance (log)", xlabelsize = 50, ylabelsize = 50, yaxisposition = :right, yticklabelalign = (:left, :center), xticklabelsvisible = false, xlabelvisible = false)
lines!(ax1, Temp_rich, everything.richness, color = ("#6B8EDE",0.8), linewidth = 5, label = "Richness")
band!(ax1, Temp_rich, everything.richness .- everything.richness_err, everything.richness .+ everything.richness_err, color = ("#6B8EDE", 0.2))
lines!(ax2, Temp_rich, everything.ϵ_var_mean, color = ("#EF8F8C", 0.8), linewidth = 5, label = "CUE Variance")
band!(ax2, Temp_rich,  everything.ϵ_var_mean .- everything.ϵ_var_err, everything.ϵ_var_mean .+ everything.ϵ_var_err, color = ("#EF8F8C", 0.2))
linkxaxes!(ax1,ax2)
l1 = [LineElement(color = ("#6B8EDE",0.8), linestyle = nothing, linewidth = 5)]
l2 = [LineElement(color = ("#EF8F8C", 0.8), linestyle = nothing, linewidth = 5)]
Legend(f[1,1], [l1, l2], tellheight = false, tellwidth = false, ["Richness", "CUE Variance"], halign = :left, valign = :top)
f

save("../result/11.png", f) 

# save("../result/rich_temp.png", f) 

# f = Figure(fontsize = 32, resolution = (1200, 900));
# ax = Axis(f[1,1], xlabel = "Temperature (°C)", ylabel = "CUE")
# lines!(ax, Temp_rich, everything.ϵ_sur_mean, color = ("#EF8F8C",1), linewidth = 5, label = "Survivor")
# band!(ax, Temp_rich, everything.ϵ_sur_mean .- everything.ϵ_sur_err , everything.ϵ_sur_mean .+ everything.ϵ_sur_err , color = ("#EF8F8C", 0.2))
# lines!(ax, Temp_rich, everything.ϵ_ext_mean, color = ("#4F363E", 0.6), linewidth = 5, label = "Extinct")
# band!(ax, Temp_rich,  everything.ϵ_ext_mean .- everything.ϵ_ext_err, everything.ϵ_ext_mean .+ everything.ϵ_ext_err, color = ("#4F363E", 0.2))
# axislegend(position = :rb)
# f
# save("../result/CUE_temp.png", f) 
