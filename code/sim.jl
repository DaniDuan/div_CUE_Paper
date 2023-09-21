# cd("WORKING_DIRECTORY")

# Load libraries
using Distributions
using LinearAlgebra
using DifferentialEquations
using Plots, StatsPlots
using Sundials
using Parameters
using CSV, DataFrames
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
richness = Float64[]; richness_var = Float64[]
for i in range(0, stop = 25, length = 26)
    T = 273.15 + i 
    # all = Float64[]
    for j in 1:50
        p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=0.4, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
        x0 = fill(0.1, (N+M))
        prob = ODEProblem(dxx!, x0, tspan, p)
        sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
        bm = sol.u[length(sol.t)][1:N]
        push!(all, length(bm[bm.>1e-7]))
        print(i, " °C Completed ",j*2, "% \n")
    end
    # rich = mean(all);rich_var = var(all)
    # push!(richness, rich); push!(richness_var, rich_var)
    # print(i, " °C Complete, ", "richness ",rich,"\n")
end 

# richness
Temp = Any[]
[append!(Temp, fill(i, 50)) for i in 0:25][1]
plot_scatter =scatter(Temp, all, color = :forestgreen, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.7)
# richness_all = DataFrame(richness = all)
# CSV.write("../data/all.csv", richness_all, writeheader=false)
# richness_all = CSV.read("../data/all.csv", DataFrame, header=false)[:,1]

################## reading EMP data ###############################
EMP_data = CSV.read("../data/EMP_filtered.csv", DataFrame, header= true)
EMP_relative_rich = EMP_data[:, "relative_rich"]
EMP_temp = EMP_data[:, "Temp"]

EMP = [EMP_temp EMP_relative_rich]
scatter!(EMP_temp, EMP_relative_rich, color = :darkorange, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.7, label = "EMP")

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
    quantile_99_richness = quantile(bin_data.Richness, 0.99)
    # Filter rows with richness above 99% quantile
    filtered_rows = filter(row -> row.Richness > quantile_99_richness, bin_data)
    append!(filtered_df, filtered_rows)
end
println(filtered_df)
################################## fitting EMP to gaussian #################################
using(LsqFit)

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
    quantile_99_richness = quantile(bin_data.Richness, 0.99)
    # Filter rows with rich above 99% quantile
    filtered_rows = filter(row -> row.Richness > quantile_99_richness, bin_data)
    append!(filtered_df_pre, filtered_rows)
end

x_data_pre = filtered_df_pre[:, "Temp"]; y_data_pre = filtered_df_pre[:, "relative_richness"]
# fitting gaussian
fit_result = curve_fit(gaussian, x_data_pre, y_data_pre, initial_params)

best_fit_params = coef(fit_result)
std_errors = stderror(fit_result)
# Generate x y values for plot
x_fit_pre = range(minimum(x_data_pre), stop=maximum(x_data_pre), length=100)
y_fit_pre = gaussian(x_fit_pre, best_fit_params)


########################################## all plots #####################
scatter(Temp, relative_all, color = :green, markerstrokecolor = nothing, 
    markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.5,label = "Prediction")
scatter!(EMP_temp, EMP_relative_rich, color = :orange, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)),  markersize = 5, fillalpha = 0.5, label = "EMP")
scatter!(x_data, y_data, color = :darkorange, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)), fillalpha = 0.7, label="")
# Plot the fitted Gaussian curve
plot!(x_fit, y_fit, linewidth=3,color = :darkorange, label="EMP_Gaussian")
scatter!(x_data_pre, y_data_pre, color = :darkgreen, markerstrokecolor = nothing, markershape = Shape(Plots.partialcircle(0, 2π)), fillalpha = 0.7, label="")
# Plot the fitted Gaussian curve
plot!(x_fit_pre, y_fit_pre, linewidth=3,color = :darkgreen, label="MiCRM_Gaussian")
xlabel!("Temperature")
ylabel!("Relative richness")

# savefig("test.pdf")

########################################### plot against EMP ##################################################
