"""
    Flux!(flu, con, p)

Given the fluid variables `flu` computes the corresponding fluxes.

As parameters (p) have: 
    χ = Function parameters

# Examples
```julia-repl
julia> 
χ = [- 1.; - 2.; - 10.]
flu = [- 5.; 0.5; 2.1; 0.5; 5.1]
Flu = [5.921317694643806; 6.02302807825841; 3.713391378258412;  4.136735467078638; 3.444079555898864]

Flux(flu,χ) - Flu # should vanish
```
"""
function Flux(flu,χ)
    μ = -flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (μ)^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    γ = (1 - v^2)^(-1//2)
    τ = -2χ₁ * x3 * T / (γ*μ^3) - 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7/5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2 - 6χ₁*x1/(γ * μ^4 * T) - 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(v^2-1)^2*x3^2)/(μ^5 * γ^2)
    Q = -10χ₀ * x2 * T / μ^3 - 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)

    f₁ = 4//3 * ρ*γ^2*v + γ*Q*(1+v^2)+ τ*v
    f₂ = 4//3 * ρ*(γ^2*v^2 + 1//4) + 2v*γ*Q + τ
    f₃ = -χ₁*γ*v*(6*γ^2 - 1)/μ^3/T - 12χ₂*(v*(6γ^2 - 1)*x1 + (6γ^2*(2v^2 + 1)-1)*x2 + v*(v^2+2)*x3)/μ^4
    f₄ = -χ₁*γ*(6γ^2*v^2+1)/μ^3/T - 12χ₂*((6γ^2*v+1)*x1 + v*(6γ^2*(1+v^2)+1)*x2 + (2v^2+1)*x3)/μ^4
    f₅ = -3χ₁*γ*v*(2γ^2*v^2+1)/μ^3/T - 12χ₂*(v*(6γ^2*v^2+3)*x1 +3*(1+6γ^2*v^2)*x2 + 3v*x3)/μ^4
    return [f₁;f₂;f₃;f₄;f₅]
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
function Is!(flu,Is,par)
    λ, κ = par
    μ = -flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (μ)^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    γ = (1 - v^2)^(-1//2)
    Is[1] = -2//5*(γ^2-1//4)*T*x1/(γ*λ) - 2γ*v*x2/T/κ - v^2*T*x3/(γ*λ)
    Is[2] = -2//5*γ*v*T*x1/λ - γ*(v^2+1)*x2/T/κ - v*T*x3/(γ*λ)
    Is[3] = -2//5*(γ^2*v^2+1//4)*T*x1/λ/γ - 2γ*v*x2/T/κ - T*x3/(γ*λ)
    return Is
end