#import Pkg; 
#Pkg.activate("Conf_Fluids")
#Pkg.add("Symbolics")
using Symbolics
using StaticArrays

function NR_step!(F, Jac, r, u0, y, p)
#    u = u0 - Jac(u0,y, p) \ F(u0,y,p)
    u0 = u0 - Jac(r,u0,y,p) \ F(r,u0,y,p)
end

"""
    c_to_f!(u, p)

Compute the fluid variables `flu` from the conservative `con` ones.
It inverts the function `F(r,flu,con,χ) = r = con(flu,χ) - con` using 
Newto-Raphson method.

As parameters have: 
    χ = Function parameters, 
    tol = error tolerance for the inversion, 
    iter_max = maximum iteration number for the NR method, 
    N = number of variables (5 here), 
    M = dimension of vector of variables

# Examples
```julia-repl
julia> 
tol = 10^(-16)
iter_max = 40
χ = [-1.0; -1.0; -5.0]
M=1
N=10
p = (χ, tol, iter_max, N, M)

u = zeros(N*M)
u0 = zeros(N*M)

con = view(reshape(u,(M,N)),:,1:5)
flu = view(reshape(u,(M,N)),:,6:10)

@time c_to_f!(u0, p) - u


```
"""
function c_to_f!(u, p)
    χ, tol, iter_max, N, M, F, Jac = p 
    #reshape(u,(N,U))
    tol = tol^2
    con = view(reshape(u,(M,N)),:,1:N÷2)
    flu = view(reshape(u,(M,N)),:,N÷2+1:N)
    for j ∈ 1:M
        r = ones(N÷2)
        iter = 1
#        while F(flu[j,:],con[j,:], χ)'*F(flu[j,:],con[j,:], χ) > tol && iter < iter_max # para F que no cambia valores
        while r'*r > tol && iter < iter_max
            flu[j,:] = NR_step!(F, Jac, r, flu[j,:], con[j,:], χ)
#            println(r)
            iter = iter + 1
            if iter == iter_max 
            println("iter_max reached j = $j")
            end
        end
    end
    return u[:]
end



"""
    f_to_c!([con, flu], p)

Compute the conservative `con` variables for the fluid ones `flu`.
It uses the function `F(flu,con,χ) = con(flu,χ) - con` with `con = 0`.

As parameters have: 
    χ = Function parameters
    N = number of variables (5 here), 
    M = dimension of vector of variables

# Examples
```julia-repl
julia> 

```
"""
function f_to_c!(u, p)
    χ, N, M, F= p
    con = view(reshape(u,(M,N)),:,1:N÷2)
    flu = view(reshape(u,(M,N)),:,N÷2+1:N)
    y = zeros(N÷2)
    yy = zeros(N÷2)
    for j ∈ 1:M
        con[j,:] = F(yy,flu[j,:],y, χ) # versión vieja
        #F(con[j,:],flu[j,:],y, χ)
    end
    return u[:]
end

function c_to_f_direct!(u,p)
    N, M = p
    c = view(reshape(u,(M,N)),:,1:N÷2)
    f = view(reshape(u,(M,N)),:,N÷2+1:N)
    
    for j ∈ 1:M
    
        e = c[j,1] # 4//3* ρ*(γ^2 - 1//4)
        s = c[j,2] # 4//3* ρ*γ^2*v
        x = s/e
        v = 3//2 * x / (1 + sqrt(1 - 3//4 * x^2))
        γ = (1 - v^2)^(-1//2)
        #e = 4//3* ρ*(γ^2 - 1/4)
        #ρ = -6 / μ^2
        ρ = 3//4 * e / (γ^2 - 1//4)
        μ = -sqrt(6/ρ)
        f[j,1] = μ    # esto es -μ
        f[j,2] = v
    end
    return u[:]
end




function f2c!(R,flu,con,χ) 
    T = (abs(flu[1]))^(-1/2) 
    v = flu[2]
    ν = flu[3]
    r1 = flu[4]
    t11 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    γ = (1 - v^2)^(-1//2)
    T00_0 = -8T^4*χ₀*(γ^2 - 1//4) #signo cambiado
    T01_0 = -8T^4*χ₀*γ^2*v
    T00_1 = χ₁*T^6*(-80//3*T^2*ν*(γ^2 - 1//4) - 20*T*r1*γ*v        - 2*t11*v^2) #signo cambiado
    T01_1 = χ₁*T^6*(-80//3*T^2*ν*γ^2*v          - 10*T*r1*γ*(1+v^2)  - 2*t11*v  )

    T00_2 = χ₂*(T^12*ν^2*(832//3+γ^2*(-1296+544*γ^2)-256//3*γ^6*(1+v^4)+γ^2*v^2*(304//3-1664//3*γ^2+512//3*γ^4+32//3*γ^2*v^2)) + T^11*r1*ν*γ*v*(-1120) + T^10*(r1^2*(96*(1+γ^4+γ^4*v^4)-288*v^2-612*γ^2+648*γ^2*v^2-192*γ^4*v^2-36*γ^2*v^4)-112*v^2*t11*ν) + T^9*r1*t11*v*γ*(v^2-1)*336 + T^8*t11^2*(9+v^2*(-42+33*v^2)+γ^2*(144*v^2-72*(1+v^4)))) #versión Marcelo

    T01_2 = χ₂*(T^12*ν^2*γ^2*v*(-3584//3 + (1600//3)γ^2*(1-v^2) - (256//3)γ^4*(1+v^4) + (512//3)*γ^4*v^2) + T^11*(-720*r1*γ*ν*(1+v^2) + 160*γ^3*r1*ν*(1-v^4)) + T^10*v*(-192*r1^2*(1+γ^4*v^2) - 112*t11*ν - 576*γ^2*r1^2*(1-v^2) + 96*γ^4*r1^2*(1+v^4)) + T^9*168*r1*t11*γ*(v^4-1) + T^8*t11^2*v*(24*(v^2-1) - 72*γ^2*(1+v^4) + 144*γ^2*v^2)) #versión Marcelo
    
    A000_1 = 3χ₁*T^5*γ*(1 - 2γ^2) #OK
    A001_1 = χ₁*T^5*γ*v*(1 - 6γ^2)  #OK 
    A011_1 = -χ₁*T^5*γ*(1 + 6γ^2*v^2) #OK
    A000_2 = -12χ₂*T^8*(T*10ν*γ*  (2γ^2-1)        + 3r1*v*(6γ^2-1)              + t11*(3v^2*γ/T)   )
    A001_2 = -12χ₂*T^8*(T*10ν*γ*v*(2γ^2-1//3)     + r1*(12γ^2*v^2  + 6γ^2 - 1)  + t11*v*γ/T*(v^2+2))
    A011_2 = -12χ₂*T^8*(T*10ν*γ*  (2γ^2*v^2+1//3) + r1*v*(6γ^2*v^2 + 12γ^2 + 1) + t11*γ/T*(2v^2+1) ) 
    
    R[1] = - con[1] + T00_0 + T00_1 + T00_2 #T00
    R[2] = - con[2] + T01_0 + T01_1 + T01_2 #T01
    R[3] = - con[3] + A000_1  + A000_2  #A000
    R[4] = - con[4] + A001_1  + A001_2  #A001
    R[5] = - con[5] + A011_1  + A011_2 #NA2_011 #A011_2  #A011
    return R[:]
end

@variables f[1:5], c[1:5], p[1:3], r[1:5]

JS_alt = Symbolics.jacobian(f2c!(r,f,c,p),f);
J_exp_alt = Symbolics.build_function(JS_alt, r, f, c, p);
Jf2c = eval(J_exp_alt[1]);


function flu_to_abs(flu)
    # cambiamos el signo de \xi_{1} y $\xi_{01}
    #s = -flu[1]
    #T = 1/sqrt(s)
    T = (abs(flu[1]))^(-1/2) 
    v = flu[2]
    ν = flu[3]
    r1 = flu[4]
    t11 = flu[5]
    γ = (1 - v^2)^(-1//2)
    #X_00 = ν*T^2/3*(4γ^2 - 1) -2r1*v*γ*T + t11*v^2
    #X_01 = 4//3 * γ^2 * v * ν*T^2 - r1*γ*T *(1 + v^2) + t11*v # Indices arriba
    #X_11 = (4 * γ^2*v^2 + 1)//3 * ν*T^2 - 2r1*v*γ*T + t11
    return [-γ/T; γ*v/T; 1.0*(ν*T^2*(4γ^2 - 1)/3-2r1*v*γ*T+t11*v^2); -1.0*(4/3*γ^2*v*ν*T^2-r1*γ*T*(1 + v^2)+t11*v); 1.0*((4*γ^2*v^2 + 1)/3*ν*T^2-2r1*v*γ*T+t11)]
end

Jf2a = Symbolics.jacobian(flu_to_abs(f),f);
Jf2a_exp = Symbolics.build_function(Jf2a, f);
Jf2a_exp1 = eval(Jf2a_exp[1]);

S = [1 0 0 0 0 ;
     0 1 0 0 0;
     0 0 3//4 0 1//4;
     0 0 0 1//2 0;
     0 0 1//4 0 3//4]