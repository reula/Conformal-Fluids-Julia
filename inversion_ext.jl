#import Pkg; 
#Pkg.activate("Conf_Fluids")
#Pkg.add("Symbolics")
using Symbolics
using StaticArrays

function F(flu,con,χ)    
    μ = flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (abs(μ))^(-1/2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    x1 = flu[3]
    x2 = flu[4]
    x3 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    
    γ = (1 - v^2)^(-1//2)
    τ =                        2χ₁ * x3 * T / (γ*μ^3)     + 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7//5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2           - 6χ₁*x1/(γ * μ^4 * T)       + 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(1-v^2)*x3^2)/(μ^5 * γ^2)
#    Q =                       10χ₁ * x2 * T / μ^3         + 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)
    Q =                       10χ₁ * x2 * T / μ^3         + 168χ₂ * x2 * (x1 + x3/γ)/(γ * μ^5)

    A = (-12 * χ₂/μ^4)*[3(2γ^2 - 1)      3v*(6γ^2 - 1)           3v^2 ;
                       v*(6γ^2 - 1)      (6γ^2*(1 + 2v^2) - 1)   v*(v^2 + 2) ;
                       (6γ^2 * v^2 + 1)  v*(6γ^2*(2 + v^2) - 1)  (2v^2 + 1)]

    e = 4//3*ρ*(γ^2 - 1//4) + 2v*γ*Q         + v^2*τ
    s = 4//3*ρ*γ^2*v        + (v^2+1//1)γ*Q + v*τ;

    Y1 = χ₁*γ/T/μ^3 *(6γ^2 - 3//1)
    #Y2 = 3χ₁*γ/T/μ^3 *v*(6γ^2 - 1//1)
    Y2 = χ₁*γ/T/μ^3 *v*(6γ^2 - 1//1)
    #Y3 = 3χ₁*γ/T/μ^3 *(6γ^2*v^2 + 1//1)
    Y3 = χ₁*γ/T/μ^3 *(6γ^2*v^2 + 1//1)
    
    return  [- con[1] + 4//3* ρ*(γ^2 - 1//4) + 2//1* Q*γ*v + τ*v^2;
             - con[2] + 4//3* ρ*γ^2*v + Q*γ*(v^2 + 1//1) + τ*v;
             - con[3] + A[1,1]*x1 + A[1,2]*x2 + A[1,3]*x3 + Y1;
             - con[4] + A[2,1]*x1 + A[2,2]*x2 + A[2,3]*x3 + Y2;
             - con[5] + A[3,1]*x1 + A[3,2]*x2 + A[3,3]*x3 + Y3]
end


@variables f[1:5], c[1:5], p[1:3]

JS = Symbolics.jacobian(F(f,c,p),f);
J_exp = Symbolics.build_function(JS, f, c, p);
Jac = eval(J_exp[1]);

function NR_step!(F, Jac, u0, y, p)
    u = u0 - Jac(u0,y, p) \ F(u0,y,p)
end

"""
    c_to_f!(u, p)

Compute the fluid variables `flu` from the conservative `con` ones.
It inverts the function `F(flu,con,χ) = con(flu,χ) - con` using 
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

flu[1,:] = [-10.0, 0.2, 1.1, 1.5, 1.2]
con[1,:] = [0.1366314976448222, 0.07009306769467444, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818]

con0 = view(reshape(u0,(M,N)),:,1:5)
flu0 = view(reshape(u0,(M,N)),:,6:10)

flu0[1,:] = [-10.0, 0.2, 1.1, 1.5, 1.2] + 0.01 .*rand(N÷2)
con0[1,:] = [0.1366314976448222, 0.07009306769467444, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818]

@time c_to_f!(u0, p) - u

10-element Array{Float64,1}:
  0.0
  0.0
  0.0
  0.0
  0.0
  4.806112841038157e-9
  2.941994980965035e-10
 -2.7242426092755068e-9
 -2.2971502477986405e-9
 -1.4237011569662172e-10
```
"""
function c_to_f!(u, p)
    χ, tol, iter_max, N, M = p 
    #reshape(u,(N,U))
    tol = tol^2
    con = view(reshape(u,(M,N)),:,1:N÷2)
    flu = view(reshape(u,(M,N)),:,N÷2+1:N)
    for j ∈ 1:M
        #flu[j,1] = -flu[j,1]
        iter = 1
        while F(flu[j,:],con[j,:], χ)'*F(flu[j,:],con[j,:], χ) > tol && iter < iter_max
            flu[j,:] = NR_step!(F, Jac, flu[j,:], con[j,:], χ)
            iter = iter + 1
            if iter == iter_max 
            println("iter_max reached j = $j")
            end
        end
    #println(iter)
        #flu[j,1] = -flu[j,1]
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
N=10
M=1
χ = [-1.0; -1.0; -5.0]
p = (χ, N, M)
u = zeros(N*M)
u0 = zeros(N*M)

con = view(reshape(u,(M,N)),:,1:5)
flu = view(reshape(u,(M,N)),:,6:10)

flu[1,:] = [-10.0, 0.2, 1.1, 1.5, 1.2]
con[1,:] = [0.1366314976448222, 0.07009306769467444, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818]

con0 = view(reshape(u0,(M,N)),:,1:5)
flu0 = view(reshape(u0,(M,N)),:,6:10)

flu0[1,:] = [-10.0, 0.2, 1.1, 1.5, 1.2]
#con0[1,:] = [0.1366314976448222, 0.07009306769467444, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818]

println(u0)
println(u)
println(f_to_c!(u0,p) - u)
println(u0)

[0.0, 0.0, 0.0, 0.0, 0.0, -10.0, 0.2, 1.1, 1.5, 1.2]
[0.1366314976448222, 0.07009306769467444, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818, -10.0, 0.2, 1.1, 1.5, 1.2]
[0.0, -1.3877787807814457e-17, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
[0.1366314976448222, 0.07009306769467442, 0.06115332989597844, 0.07178418128379448, 0.04927907295689818, -10.0, 0.2, 1.1, 1.5, 1.2]
```
"""
function f_to_c!(u, p)
    χ, N, M = p
    con = view(reshape(u,(M,N)),:,1:N÷2)
    flu = view(reshape(u,(M,N)),:,N÷2+1:N)
    y = zeros(N÷2)
    for j ∈ 1:M
        #flu[j,1] = -flu[j,1]
        con[j,:] = F(flu[j,:],y, χ)
        #flu[j,1] = -flu[j,1]
    end
    return u[:]
end

function c_to_f_direct!(u,p)
    χ, N, M = p
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

"""
    split_u!(u,c,f,N,M) NO FUNCIONA!

split in two matrices the vector u according to the sizes given.
This is a view, not a copy!

    N = number of variables (5x2 here), 
    M = dimension of vector of variables

# Examples
```julia-repl
julia> 
N=10
M=1

u = zeros(N*M)

```
"""
function split_u!(u,c,f,N,M)
    c = view(reshape(u,(M,N)),:,1:N÷2)
    f = view(reshape(u,(M,N)),:,N÷2+1:N)
end



function F_alt(flu,con,χ)    
    μ = flu[1]  # esto es -μ
    #μ = view(flu,1)
    T = (abs(μ))^(-1/2) # use μ positive, so I changed μ -> -μ
    v = flu[2]
    ν = flu[3]
    r1 = flu[4]
    t11 = flu[5]
    χ₀ = χ[1]
    χ₁ = χ[2]
    χ₂ = χ[3]
    
    γ = (1 - v^2)^(-1//2)
    τ =                        2χ₁ * x3 * T / (γ*μ^3)     + 24χ₂*(1//2*(1-v^2)x3^2 + 14//3 * x2^2 + 7//5*(1-v^2)x1*x3)/μ^5
    ρ = -6χ₀ / μ^2           - 6χ₁*x1/(γ * μ^4 * T)       + 42χ₂*(6//5 *x1^2 + 10γ^2*x2^2 + 3//2*(1-v^2)*x3^2)/(μ^5 * γ^2)
#    Q =                       10χ₁ * x2 * T / μ^3         + 168χ₂ * x2 * (x1 - (v^2 - 1)x3)/(γ * μ^5)
    Q =                       10χ₁ * x2 * T / μ^3         + 168χ₂ * x2 * (x1 + x3/γ)/(γ * μ^5)

    e = 4//3*ρ*(γ^2 - 1//4) + 2v*γ*Q         + v^2*τ
    s = 4//3*ρ*γ^2*v        + (v^2+1//1)γ*Q + v*τ;
    
    return  [- con[1] + 4//3* ρ*(γ^2 - 1//4) + 2//1* Q*γ*v + τ*v^2;
             - con[2] + 4//3* ρ*γ^2*v + Q*γ*(v^2 + 1//1) + τ*v;
             - con[3] -3χ₁*T^5*γ*(2γ^2-1)   - 12χ₂*(T*10ν*γ*  (2γ^2-1)      
                                                    + 3r1*v(6γ^2-1)
                                                    + t11*(3v^2*γ/T)); #A000
             - con[4] -χ₁*T^5*γ*v(6γ^2-1)   - 12χ₂*(T*10ν*γ*v*(2γ^2-1//3)
                                                    + r1*(12γ^2*v^2 + 6γ^2 - 1)
                                                    + t11*v*γ/T*(v^2+2)); #A001
             - con[5] -χ₁*T^5*γ(6γ^2*v^2+1) - 12χ₂*(T*10ν*γ*  (2γ^2*v^2+1//3)
                                                    + r1*v*(6γ^2*v^2 + 12γ^2 + 1)
                                                    + t11*γ/T*(2v^2+1)] #A011
end

