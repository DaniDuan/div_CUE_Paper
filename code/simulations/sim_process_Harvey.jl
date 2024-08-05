include("./sim_frame.jl")

using Glob
path = glob("iters0_*", "../data/output_0/")

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

everything = zeros(Float64, num_temps, 44)

@time for j in 1: num_temps
    all_rich_H = Float64[]; ϵ_sur_H = Float64[]; ϵ_ext_H = Float64[]; all_ϵ_H = Float64[]; all_ϵ_r0_H = Float64[]; 
    feas_H = Float64[]; ELV_H = Float64[]; all_Shannon_H =  Float64[]; all_Simpson_H = Float64[]; #Eϵ_sur = Float64[]; Eϵ_ext = Float64[]; 
    ϵ_var_H = Float64[]; u_sur_H = Float64[]; u_ext_H = Float64[]; m_sur_H = Float64[]; m_ext_H = Float64[];
    all_Eu_H = Float64[]; all_Em_H = Float64[]; all_Eu_sur_H = Float64[]; all_Em_sur_H = Float64[]; 
    all_Tpu_H = Float64[]; all_Tpm_H = Float64[]; all_Tpu_sur_H = Float64[]; all_Tpm_sur_H = Float64[]
    
    for i in 1:length(path)
        @load path[i] all_rich ϵ_sur ϵ_ext all_ϵ all_ϵ_r0 feas ELV all_Shannon all_Simpson ϵ_var u_sur u_ext m_sur m_ext all_Eu all_Em all_Eu_sur all_Em_sur all_Tpu all_Tpm all_Tpu_sur all_Tpm_sur
        
        push!(all_rich_H, all_rich[j]); push!(feas_H, feas[j]); push!(ELV_H, ELV[j]);
        push!(all_Shannon_H, all_Shannon[j]); push!(all_Simpson_H, all_Simpson[j]);
        append!(ϵ_sur_H, ϵ_sur[j]); append!(ϵ_ext_H, ϵ_ext[j]); 
        append!(all_ϵ_H, all_ϵ[j]); push!(ϵ_var_H, ϵ_var[j]); append!(all_ϵ_r0_H, all_ϵ_r0[j]); #CUE
        append!(all_Eu_H, all_Eu[j]); append!(all_Em_H, all_Em[j]); append!(all_Eu_sur_H, all_Eu_sur[j]); append!(all_Em_sur_H, all_Em_sur[j]); #Eu Em
        append!(all_Tpu_H, all_Tpu[j]); append!(all_Tpm_H, all_Tpm[j]); append!(all_Tpu_sur_H, all_Tpu_sur[j]); append!(all_Tpm_sur_H, all_Tpm_sur[j]); #Tp
        append!(u_sur_H, u_sur[j]); append!(u_ext_H, u_ext[j]);
        append!(m_sur_H, m_sur[j]); append!(m_ext_H, m_ext[j])
    end 
    everything[Int(j),:] = [mean(all_rich_H), std(all_rich_H)/sqrt(length(all_rich_H)), mean(feas_H), std(feas_H)/sqrt(length(feas_H)), mean(ELV_H), std(ELV_H)/sqrt(length(ELV_H)),
    mean(all_Shannon_H), std(all_Shannon_H)/sqrt(length(all_Shannon_H)), mean(all_Simpson_H), std(all_Simpson_H)/sqrt(length(all_Simpson_H)),
    mean(ϵ_sur_H), std(ϵ_sur_H)/sqrt(length(ϵ_sur_H)), mean(ϵ_ext_H), std(ϵ_ext_H)/sqrt(length(ϵ_ext_H)), 
    mean(all_ϵ_H), std(all_ϵ_H)/sqrt(length(all_ϵ_H)), mean(ϵ_var_H), std(ϵ_var_H)/sqrt(length(ϵ_var_H)), mean(all_ϵ_r0_H), std(all_ϵ_r0_H)/sqrt(length(all_ϵ_r0_H)),
    mean(all_Eu_H), std(all_Eu_H)/sqrt(length(all_Eu_H)), mean(all_Em_H), std(all_Em_H)/sqrt(length(all_Em_H)),
    mean(all_Eu_sur_H), std(all_Eu_sur_H)/sqrt(length(all_Eu_sur_H)), mean(all_Em_sur_H), std(all_Em_sur_H)/sqrt(length(all_Em_sur_H)),
    mean(all_Tpu_H), std(all_Tpu_H)/sqrt(length(all_Tpu_H)), mean(all_Tpm_H), std(all_Tpm_H)/sqrt(length(all_Tpm_H)),
    mean(all_Tpu_sur_H), std(all_Tpu_sur_H)/sqrt(length(all_Tpu_sur_H)), mean(all_Tpm_sur_H), std(all_Tpm_sur_H)/sqrt(length(all_Tpm_sur_H)),
    mean(u_sur_H), std(u_sur_H)/sqrt(length(u_sur_H)), mean(u_ext_H), std(u_ext_H)/sqrt(length(u_ext_H)), 
    mean(m_sur_H), std(m_sur_H)/sqrt(length(m_sur_H)), mean(m_ext_H), std(m_ext_H)/sqrt(length(m_ext_H))]
    print(j-1, " °C Complete, ", "richness ", mean(all_rich_H),"\n") 

end

col_names = ["richness", "richness_err", "feas", "feas_err", "ELV", "ELV_err", "Shannon", "Shannon_err", "Simpson", "Simpson_err",
            "ϵ_sur_mean", "ϵ_sur_err", "ϵ_ext_mean", "ϵ_ext_err", "ϵ", "ϵ_err","ϵ_var_mean", "ϵ_var_err", "ϵ_r0", "ϵ_r0_err",
            "Eu", "Eu_err", "Em", "Em_err", "Eu_sur", "Eu_sur_err", "Em_sur", "Em_sur_err",
            "Tpu", "Tpu_err", "Tpm", "Tpm_err", "Tpu_sur", "Tpu_sur_err", "Tpm_sur", "Tpm_sur_err", 
            "u_sur_mean", "u_sur_err", "u_ext_mean", "u_ext_err", "m_sur_mean", "m_sur_err", "m_ext_mean", "m_ext_err"];
everything = DataFrame(everything, col_names);

### Saving results
# CSV.write("../data/temp_gradient.csv", everything, writeheader=false)
