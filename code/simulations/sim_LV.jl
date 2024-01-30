include("../sim_frame.jl");

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

T = 273.15+15
p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=L, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
## run simulation
prob = ODEProblem(dxx!, x0, tspan, p)
sol =solve(prob, AutoVern7(Rodas5()), save_everystep = false, callback=cb)

@time p_lv = Eff_LV_params(p=p, sol=sol);
p_lv.ℵ
# # Set initial conditions
# Ci = fill(0.0, N)
# for i in 1:N
#     Ci[i] = 0.1
# end

# # Define ODE problem
# # prob_LV = ODEProblem(LV_dx!, Ci, (0.0, t_span), LV1)
# t_span = 15000
# prob_LV = ODEProblem(LV_dx!, Ci, (0.0, t_span), p_lv)

# # Integrate ODE with callback and terminate if it fails to converge more than 3 times (this is to avoid unncesessary overhead)
# sol_LV = solve(prob_LV, reltol=1e-9,abstol=1e-9, CVODE_BDF(max_convergence_failures = 3, max_hnil_warns = 1, stability_limit_detect=true), callback=cb)

# # Calculate effective LV jacobian
# LV_jac = Eff_Lv_Jac(p_lv=p_lv, sol=sol)
