"""
    Flux(Fl,u, p)

Given the fluid variables `flu` and `con` computes the corresponding fluxes.
The value of `con` is used to approximate the next value of it when computed from `flu`. 

As parameters (p) have: 
    χ = Function parameters


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

function Flux_imp!(Fl,u,p)
    flu = u[6:end]
    con = u[1:5]
    χ = p
    T = (abs(flu[1]))^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    ν = flu[3]
    r1 = flu[4]
    t11 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    γ = (1 - v^2)^(-1//2)
    T11_0 = -8χ₀*T^4*(γ^2*v^2+1//4)
    T11_1 = χ₁*(-80//3*T^8*ν*(γ^2*v^2+1//4)  - 20*T^7*r1*γ*v  - 2*T^6*t11)
    T11_2 = χ₂*T^8*(ν^2*(-560*T^4*(γ^2*v^2 + 1//4)) - ν*r1*T^3*1120*γ*v - ν*t11*T^2*112
                    + r1^2*T^2*(-315*v^2-217)
                    + r1*t11*T*336*γ*v*(v^2-1)
                    + t11^2*(-111 - 78v^2 + 189v^4)/4)
    #A111_1 = -3T^5*γ*χ₁*v*(1 + 2γ^2*v^2) OK
    Fl[1] = con[2]
    Fl[2] = T11_0 + T11_1 + T11_2
    Fl[3] = con[4]
    Fl[4] = con[5]
    Fl[5] = -3χ₁*γ*v*T^5*(2γ^2*v^2+1) - 12χ₂*(v*T*10*γ*(2γ^2*v^2+1)*ν  + 3*(1+6γ^2*v^2)*r1  + 3v*γ*t11/T)
    Fl[6] = 0.
    Fl[7] = 0.
    Fl[8] = 0.
    Fl[9] = 0.
    Fl[10] = 0.
    return -Fl[:]
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

function Is!(sourcevec,u,t,par)
    χ, ξ = par #only used to test otherwise are defined below. a is and amplitude to change I from outside
    Is = zeros(10)
    flu = u[6:end]
    μ = flu[1] 
    T = (abs(μ))^(-1//2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    γ = (1 - v^2)^(-1//2)
    x1 = γ*T*flu[3]
    x2 = flu[4]
    x3 = γ/T*flu[5]
    sourcevec[1] = 0.0
    sourcevec[2] = 0.0
    sourcevec[3] = (-(γ^2    -1//4)*x1*ξ[1]        + 2γ*v*x2*ξ[2]/T      - μ*v^2*x3*ξ[3])/μ^5;
    sourcevec[4] = (-(γ^2*v  + 0  )*x1*ξ[1]        + γ*(v^2+1)*x2*ξ[2]/T - μ*v*x3*ξ[3])/μ^5;
    sourcevec[5] = (-(γ^2*v^2+1//4)*x1*ξ[1]        + 2γ*v*x2*ξ[2]/T      - μ*x3*ξ[3])/μ^5;
    sourcevec[6] = 0.0
    sourcevec[7] = 0.0
    sourcevec[8] = 0.0
    sourcevec[9] = 0.0
    sourcevec[10] = 0.0
end


function Is_dummy!(sourcevec,u,t,par)
    sourcevec[1] = 0.0
    sourcevec[2] = 0.0
    sourcevec[3] = 0.0
    sourcevec[4] = 0.0
    sourcevec[5] = 0.0
    sourcevec[6] = 0.0
    sourcevec[7] = 0.0
    sourcevec[8] = 0.0
    sourcevec[9] = 0.0
    sourcevec[10] = 0.0
end




function Speed_max(u, par_flux)
    #  Here we compute the maximal propagation speed of the equation, for the cases of real eigenvalues is the spectral radious of the 
    #  Jacobian (when the roots have imaginary values I guess it is the maximal real part of the eigenvalues).
    χ = par_flux
    return 1. #para revisar.... debiera ser una función de varias cosas... 
end