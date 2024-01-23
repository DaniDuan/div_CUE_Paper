include("../sim_frame.jl")

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

################### binning by integers of temp ##############
filtered_data = copy(EMP_data)
filtered_data.Temp = round.(Int, filtered_data.Temp)
# bins created on unique integer temps
bins = unique(filtered_data.Temp)
filtered_df = DataFrame()
# EMP_meanerr = zeros(Float64, num_temps, 3) # for mean plot 
# n = 0
for unique_bin in bins
    # n = n+1
    # Filter rows for bin
    bin_data = filter(row -> row.Temp == unique_bin, filtered_data)

    # #getting the average and err 
    # relative_bin_data = bin_data.Richness/maximum(bin_data.Richness)
    # mean_bin_data = mean(relative_bin_data)
    # err_bin_data = std(relative_bin_data)/sqrt(length(relative_bin_data))

    # Calculate 99% quantile of richness for bin
    quantile_97_richness = quantile(bin_data.Richness, 0.97)
    # Filter rows with richness above 99% quantile
    filtered_rows = filter(row -> row.Richness > quantile_97_richness, bin_data)
    append!(filtered_df, filtered_rows)
    # EMP_meanerr[Int(n),:] = [unique_bin, mean(relative_bin_data), std(relative_bin_data)/sqrt(length(relative_bin_data))]
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


