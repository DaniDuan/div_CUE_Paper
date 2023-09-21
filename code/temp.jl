"
To include temperature function 
"

# function temp_trait(N, T; ρ_t = [0 0], Tr= 273.15+10, Ed = 3.5, L=0.4) 
function temp_trait(N, kw) 
    k = 0.0000862 # Boltzman constant
    @unpack T, Tr, Ed= kw
    B,E,Tp = randtemp_param(N, kw)
    temp_p = B .* exp.((-E./k) * ((1/T)-(1/Tr)))./(1 .+ (E./(Ed .- E)) .* exp.(Ed/k * (1 ./Tp .- 1/T)))
    return temp_p
end

# function randtemp_param(N; ρ_t =  [-1.0 -1.0], Tr= 273.15+13, Ed= 3.5, L=0.4)
function randtemp_param(N, kw)
    @unpack T, ρ_t, Tr, Ed, L = kw
    ρ_t[ρ_t .== 1.0] .= 1-eps()
    ρ_t[ρ_t .== -1.0] .= -1+eps()
    k = 0.0000862 # Boltzman constant
    B0 = [log((0.138/(1 - L - 0.22)) * exp((-0.82/k) * ((1/Tr)-(1/273.15)))/(1 + (0.82/(Ed - 0.82)) * exp(Ed/k * (1/308.15 - 1/Tr)))) log(0.138 *exp((-0.67/k) * ((1/Tr)-(1/273.15)))/(1 + (0.67/(Ed - 0.67)) * exp(Ed/k * (1/311.15 - 1/Tr))))]# Using CUE0 = 0.22, mean growth rate = 0.48
    B0_var = 0.17*abs.(B0); E_mean = [0.82 0.67]; E_var =  0.14*abs.(E_mean)
    cov_xy = ρ_t .* B0_var.^0.5 .* E_var .^ 0.5
    meanv = [B0 ; E_mean]
    cov_u = [B0_var[1] cov_xy[1]; cov_xy[1] E_var[1]]
    cov_m = [B0_var[2] cov_xy[2]; cov_xy[2] E_var[2]]
    allu = rand(MvNormal(meanv[:,1], cov_u), N)
    allm = rand(MvNormal(meanv[:,2], cov_m), N)
    B = [exp.(allu[1,:]) exp.(allm[1,:])]
    E = [allu[2,:] allm[2,:]]
    
    Tpu = 273.15 .+ rand(Normal(35, 5), N)
    Tpm = Tpu .+ 3
    Tp = [Tpu Tpm]
    return B,E,Tp
end 


# # Load the CSV data
# data = CSV.File("../data/summary.csv") |> DataFrame
# f_data = filter(row -> row.source == "meta", data)
# var_B0_ss = var(log.(f_data.B0_ss))/mean(log.(f_data.B0_ss))
# var_E_ss = var(f_data.E_ss)/mean(f_data.E_ss)
# cov_B0_E = cov(log.(f_data.B0_ss), f_data.E_ss)
# print(var_B0_ss)
# print(var_E_ss)
# print(cov_B0_E)
