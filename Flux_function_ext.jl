"""
    Flux(con, p)

Given the fluid variables `flu` and `con` computes the corresponding fluxes.
The value of `con` is used to approximate the next value of it when computed from `flu`. 

As parameters (p) have: 
    p = (flu, χ)
    χ = Function parameters
    con = value of `con` is used to approximate the next value of it when computed from `flu`.

# Examples
```julia-repl
julia> 
χ = [- 1.; - 2.; - 10.]
con = zeros(N)
flu = [- 5.; 0.5; 2.1; 0.5; 5.1]
Flu = [5.921317694643806; 6.02302807825841; 3.713391378258412;  4.136735467078638; 3.444079555898864]

Flux(flu,χ) - Flu # should vanish
```
"""
function Flux(u,p)
    χ = p
    f = zeros(length(u))
    #c_to_f(u,χ)    
    #con = view[1:5]
    flu = u[6:end]
    μ = flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (abs(μ))^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    γ = (1 - v^2)^(-1//2)
    τ = 2χ₁ * x3 * T / (γ*μ^3) + 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7/5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2 - 6χ₁*x1/(γ * μ^4 * T) + 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(v^2-1)^2*x3^2)/(μ^5 * γ^2)
    Q = 10χ₀ * x2 * T / μ^3 + 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)

    f[1] = 4//3 * ρ*γ^2*v + γ*Q*(1+v^2)+ τ*v
    f[2] = 4//3 * ρ*(γ^2*v^2 + 1//4) + 2v*γ*Q + τ
    f[3] = χ₁*γ*v*(6*γ^2 - 1)/μ^3/T - 12χ₂*(v*(6γ^2 - 1)*x1 + (6γ^2*(2v^2 + 1)-1)*x2 + v*(v^2+2)*x3)/μ^4
    f[4] = χ₁*γ*(6γ^2*v^2+1)/μ^3/T - 12χ₂*((6γ^2*v+1)*x1 + v*(6γ^2*(1+v^2)+1)*x2 + (2v^2+1)*x3)/μ^4
    f[5] = 3χ₁*γ*v*(2γ^2*v^2+1)/μ^3/T - 12χ₂*(v*(6γ^2*v^2+3)*x1 +3*(1+6γ^2*v^2)*x2 + 3v*x3)/μ^4
    return -f[:]
end

function Flux_cut(u,p)
    χ = p
    f = zeros(length(u))
    #c_to_f(u,χ)    
    #con = view[1:5]
    flu = u[6:end]
    μ = flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (abs(μ))^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    γ = (1 - v^2)^(-1//2)
    τ = 2χ₁ * x3 * T / (γ*μ^3) + 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7/5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2 - 6χ₁*x1/(γ * μ^4 * T) + 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(v^2-1)^2*x3^2)/(μ^5 * γ^2)
    Q = 10χ₀ * x2 * T / μ^3 + 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)

    f[1] = 0.0 #4//3 * ρ*γ^2*v + γ*Q*(1+v^2)+ τ*v
    f[2] = 0.0 #4//3 * ρ*(γ^2*v^2 + 1//4) + 2v*γ*Q + τ
    f[3] = χ₁*γ*v*(6*γ^2 - 1)/μ^3/T - 12χ₂*(v*(6γ^2 - 1)*x1 + (6γ^2*(2v^2 + 1)-1)*x2 + v*(v^2+2)*x3)/μ^4
    f[4] = χ₁*γ*(6γ^2*v^2+1)/μ^3/T - 12χ₂*((6γ^2*v+1)*x1 + v*(6γ^2*(1+v^2)+1)*x2 + (2v^2+1)*x3)/μ^4
    f[5] = 3χ₁*γ*v*(2γ^2*v^2+1)/μ^3/T - 12χ₂*(v*(6γ^2*v^2+3)*x1 +3*(1+6γ^2*v^2)*x2 + 3v*x3)/μ^4
    return f[:]
end











"""
    Is!(flu, Is, p)

Given the fluid variables `flu` computes the sources `Is`.

As parameters have: 
    p = decay constants, 

# Examples
```julia-repl
julia> 
Is = zeros(3)
par = (49.735919716217296,16.578639905405765)
I_c = [- 0.09488575328013785; - 0.12155655385334033; - 0.12140081195218293]
Is!(flu,Is,par) - I_c # should be very small
```
"""
function Is(u,t,par)
    χ, ζ = par #only used to test otherwise are defined below. a is and amplitude to change I from outside
    #κ = ζ[1] # not used
    #λ = ζ[2] # not used
    a = ζ[3] 
    Is = zeros(10)
    flu = u[6:end]
    μ = flu[1] 
    χ₀= χ[1]
    χ₁= χ[2] # lo hacemos positivo
    if χ₁ < 0
        κ = χ₀*μ^5/(15π*χ₁^2) # OK 
        λ = -χ₀*μ^4/(π*χ₁^2)   # OK
    else
        κ = 1.
        λ = 1.
    end
    T = (abs(μ))^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    γ = 1. #(1 - v^2)^(-1//2)
    Is[1] = 0.
    Is[2] = 0.
    Is[3] = -2//5*a*(γ^2-1//4)*T*x1/(γ*λ) - 2γ*v*x2/T/κ - v^2*T*x3/(γ*λ)
    Is[4] = -2//5*a*γ*v*T*x1/λ - γ*(v^2+1)*x2/T/κ - v*T*x3/(γ*λ)
    Is[5] = -2//5*a*(γ^2*v^2+1//4)*T*x1/λ/γ - 2γ*v*x2/T/κ - T*x3/(γ*λ)
    return Is[:]      #(1 - ℯ^(-5. *t))
end

function Is_dummy(u,t,par)
    Is = zeros(10)
    return Is[:]
end

function Speed_max(u, par_flux)
    #  Here we compute the maximal propagation speed of the equation, for the cases of real eigenvalues is the spectral radious of the 
    #  Jacobian (when the roots have imaginary values I guess it is the maximal real part of the eigenvalues).
    χ = par_flux
    return 1. #para revisar.... debiera ser una función de varias cosas... 
end