

def_m(N, M, kw) = ones(N)
def_ρ(N, M, kw) = ones(M)
def_ω(N, M, kw) = ones(M)
def_u(N, M, kw) = copy(rand(Distributions.Dirichlet(M, 1.0), N)')

function def_l(N, M, kw)
    l = zeros(N, M, M)
    ϕ = fill(1.0, M)
    dD = Dirichlet(ϕ[:])
    for i = 1:N
        for α = 1:M
            l[i, α, :] = rand(dD) * kw[:L]
        end
    end
    return l
end

function generate_params(N, M; f_m=def_m, f_ρ=def_ρ, f_ω=def_ω, f_u=def_u, f_l=def_l, kwargs...)
    kw = Dict{Symbol, Any}(kwargs) # for temp, include: T, ρ_t, Tr, Ed, L
    if haskey(kw, :T)
        tt = temp_trait(N, kw)#[1]
        push!(kw, :tt => tt)
    end 
    # consumers
     l = f_l(N, M, kw)
     λ = zeros(N, M)
     m = f_m(N, M, kw)
     u = f_u(N, M, kw)

     ## Calculating total leakage of consumers per resource
    for i in 1:N
        for α in 1:M
            λ[i, α] = sum(l[i, α, :])
        end
    end

    # Eϵ = temp_trait(N, kw)[2]
     # resources
     ρ = f_ρ(N, M, kw)
     ω = f_ω(N, M, kw)

     kw_nt = (; kwargs...)
    #  p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ, Eϵ)
     p_nt = (N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)

     out = Base.merge(p_nt, kw_nt)

     return out #(N=N, M=M, u=u, m=m, l=l, ρ=ρ, ω=ω, λ=λ)
 end

# p = generate_params(N, M; f_u=F_u, f_m=F_m, f_ρ=F_ρ, f_ω=F_ω, L=l_α, T=T, ρ_t=ρ_t, Tr=Tr, Ed=Ed)
# p.Eϵ