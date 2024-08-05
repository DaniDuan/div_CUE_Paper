include("./sim_frame.jl")
# using JLD2

N=100
M=50
### Temp params 
ρ_t= [0.0000 0.0000]; # realistic covariance
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

# Retrieve the environment variable as a string
index_str = ENV["SLURM_ARRAY_TASK_ID"]
# Convert the string to a numeric value (e.g., Integer)
index = parse(Int, index_str)

all_rich = Float64[]; ϵ_sur = Vector{Vector{Float64}}(); ϵ_ext = Vector{Vector{Float64}}(); all_ϵ = Vector{Vector{Float64}}(); all_ϵ_r0 = Vector{Vector{Float64}}(); 
feas = Float64[]; ELV = Float64[]; all_Shannon =  Float64[]; all_Simpson = Float64[]; #Eϵ_sur = Float64[]; Eϵ_ext = Float64[]; 
ϵ_var = Float64[]; u_sur =  Vector{Vector{Float64}}(); u_ext =  Vector{Vector{Float64}}(); m_sur =  Vector{Vector{Float64}}(); m_ext =  Vector{Vector{Float64}}();
all_Eu = Vector{Vector{Float64}}(); all_Em =  Vector{Vector{Float64}}(); all_Eu_sur = Vector{Vector{Float64}}(); all_Em_sur = Vector{Vector{Float64}}(); 
all_Tpu = Vector{Vector{Float64}}(); all_Tpm =  Vector{Vector{Float64}}(); all_Tpu_sur =  Vector{Vector{Float64}}(); all_Tpm_sur =  Vector{Vector{Float64}}()

@time for i in range(0, stop = 30, length = 31)
    T = 273.15 + i 
    p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
    ## Calc CUE
    ϵ = (p.u * x0[N+1:N+M] .* (1 .- p.L) .- p.m) ./ (p.u * x0[N+1:N+M])
    ## run simulation
    prob = ODEProblem(dxx!, x0, tspan, p)
    sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm = sol.u[length(sol.t)][1:N]
    sur = (1:N)[bm .> 1.0e-7]
    ## Shannon & Simpson
    p_v = bm[sur]./sum(bm[sur])
    Shannon = - sum(p_v .* log.(p_v))
    Simpson = 1/ sum(p_v .^2)
    ## collecting E and Tp
    Eu = p.E[:,1]; Em = p.E[:,2]
    Eu_sur = Eu[sur]; Em_sur = Em[sur]
    Tpu = p.Tp[:,1]; Tpm = p.Tp[:,2]
    Tpu_sur = Tpu[sur]; Tpm_sur = Tpm[sur]
    ###
    p_lv = Eff_LV_params(p=p, sol=sol);
    ## Making prediction based on Tom's paper
    r0 = (1:N)[p_lv.r .> 0]
    N_sur = length(r0) 
    ϵ_r0 = ϵ[r0]
    sur_ℵ = p_lv.ℵ[r0, r0]  # eliminate all species with r<0
    mean_ℵ = mean([sur_ℵ[i, j]/diag(sur_ℵ)[i] for i in 1:N_sur for j in 1:N_sur if i != j])
    sur_r = p_lv.r[r0] 
    mean_r = mean(sur_r) 
    pred = sum(sur_r .> ((N_sur - 1)* mean_ℵ)*mean_r/(1+(N_sur-1)*mean_ℵ)) # upper bound for the number of survivors
    ## running LV
    prob_LV = ODEProblem(LV_dx!, Ci, tspan, p_lv)
    sol_LV = solve(prob_LV, AutoVern7(Rodas5()), save_everystep = false, callback=cb)
    bm_LV = sol_LV.u[length(sol_LV.t)][1:N]
    pred_LV = sum(bm_LV.>1e-7)

    push!(all_rich, length(sur)); push!(feas, pred); push!(ELV, pred_LV);
    push!(all_Shannon, Shannon); push!(all_Simpson, Simpson);
    push!(ϵ_sur, ϵ[sur]); push!(ϵ_ext, ϵ[bm .< 1.0e-7]); 
    push!(all_ϵ, ϵ); push!(ϵ_var, log(var(ϵ))); push!(all_ϵ_r0, ϵ_r0); #CUE
    push!(all_Eu, Eu); push!(all_Em, Em); push!(all_Eu_sur, Eu_sur); push!(all_Em_sur, Em_sur); #Eu Em
    push!(all_Tpu, Tpu); push!(all_Tpm, Tpm); push!(all_Tpu_sur, Tpu_sur); push!(all_Tpm_sur, Tpm_sur); #Tp
    push!(u_sur, sum(p.u, dims = 2)[sur]); push!(u_ext, sum(p.u, dims = 2)[bm .< 1.0e-7]);
    push!(m_sur, p.m[sur]); push!(m_ext, p.m[bm .< 1.0e-7])

end 

@save "../data/iters0_$(index).jld2" all_rich ϵ_sur ϵ_ext all_ϵ all_ϵ_r0 feas ELV all_Shannon all_Simpson ϵ_var u_sur u_ext m_sur m_ext all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur
# @load "../data/iters_1.jld2" all_rich ϵ_sur ϵ_ext all_ϵ all_ϵ_r0 feas ELV all_Shannon all_Simpson ϵ_var u_sur u_ext m_sur m_ext all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur

    # everything[Int(i+1),:] = [mean(all_rich), std(all_rich)/sqrt(length(all_rich)), mean(feas), std(feas)/sqrt(length(feas)), mean(ELV), std(ELV)/sqrt(length(ELV)),
    #     mean(all_Shannon), std(all_Shannon)/sqrt(length(all_Shannon)), mean(all_Simpson), std(all_Simpson)/sqrt(length(all_Simpson)),
    #     mean(ϵ_sur), std(ϵ_sur)/sqrt(length(ϵ_sur)), mean(ϵ_ext), std(ϵ_ext)/sqrt(length(ϵ_ext)), 
    #     mean(all_ϵ), std(all_ϵ)/sqrt(length(all_ϵ)), mean(ϵ_var), std(ϵ_var)/sqrt(length(ϵ_var)), mean(all_ϵ_r0), std(all_ϵ_r0)/sqrt(length(all_ϵ_r0)),
    #     mean(all_Eu), std(all_Eu)/sqrt(length(all_Eu)), mean(all_Em), std(all_Em)/sqrt(length(all_Em)),
    #     mean(all_Eu_sur), std(all_Eu_sur)/sqrt(length(all_Eu_sur)), mean(all_Em_sur), std(all_Em_sur)/sqrt(length(all_Em_sur)),
    #     mean(all_Tpu), std(all_Tpu)/sqrt(length(all_Tpu)), mean(all_Tpm), std(all_Tpm)/sqrt(length(all_Tpm)),
    #     mean(all_Tpu_sur), std(all_Tpu_sur)/sqrt(length(all_Tpu_sur)), mean(all_Tpm_sur), std(all_Tpm_sur)/sqrt(length(all_Tpm_sur)),
    #     mean(u_sur), std(u_sur)/sqrt(length(u_sur)), mean(u_ext), std(u_ext)/sqrt(length(u_ext)), 
    #     mean(m_sur), std(m_sur)/sqrt(length(m_sur)), mean(m_ext), std(m_ext)/sqrt(length(m_ext))]
    # print(i, " °C Complete, ", "richness ", mean(all_rich),"\n") 

